/*
 * mandelbrot synthesizer
 *
 * Copyright (C) 2019-2019 Stefan Sauer <ensonic@users.sf.net>
 *
 * gcc -Wall -g  mandelsynth3.c -o mandelsynth3 `pkg-config gtk+-3.0 cairo gstreamer-1.0 gstreamer-app-1.0 --cflags --libs` -lm
 *
 * Usage:
 * - left-mouse down + drag: change the sprectrum
 * - '1': decrement center mode
 * - '2': increment center mode
 * - '3': decrement octave
 * - '4': increment octave
 * - alpha keys (1 oct) to play notes ('y/z') -> c, 's' -> c#, ...
 *
 * need to run
 * gsettings set org.gnome.desktop.peripherals.touchpad disable-while-typing false
 * before starting this to be able to move the mouse while typing
 */
/* TODO:
 * - build array of num_harmonics until we reach SRATE/2
 *   - this will be niter for the fractal
 * - make it polyphonic
 *   - extract setting to voice object
 *   - have a simple voice allocation
 * - support morphing between two points in the fractal
 *   - e.g. left to move start-point, right to move end-point
 *   - a) over time using an envelope from buzztrax
 *   - b) using 2 lfos for x/y
 *   - c) along some path
 * - x/y or polar coords?
 * - try initializing the phase of the osc from the fractal values
 * - try slightly detuning the harmonics
 *
 * - turn into vcvrack module to prototype interface and reuse modulators
 *   - have 2 osc per voice
 *   - each osc has {x1,x,x2},{y1,y,y2} that can be modulated
 *     - x1, x2 set the range, the x target is 0...100 and picks a spectrum in
 *      the range, likewise for y
 *   - offer various combine modes for the oscs (mix, am, fm, ...)
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>
#include <gtk/gtk.h>
#include <gst/gst.h>
#include <gst/app/gstappsrc.h>
#include <gst/audio/audio.h>
#include <gst/fft/gstfftf64.h>

#define SRATE 44100

// generated with midi2frq (octaves 1-7)
float midi2frq[]={
    32.703,    34.648,    36.708,    38.891,    41.203,    43.654,    46.249,    48.999,    51.913,    55.000,    58.270,    61.735,
    65.406,    69.296,    73.416,    77.782,    82.407,    87.307,    92.499,    97.999,   103.826,   110.000,   116.541,   123.471,
   130.813,   138.591,   146.832,   155.563,   164.814,   174.614,   184.997,   195.998,   207.652,   220.000,   233.082,   246.942,
   261.626,   277.183,   293.665,   311.127,   329.628,   349.228,   369.994,   391.995,   415.305,   440.000,   466.164,   493.883,
   523.251,   554.365,   587.329,   622.254,   659.255,   698.456,   739.989,   783.991,   830.609,   880.000,   932.328,   987.766,
  1046.502,  1108.731,  1174.659,  1244.508,  1318.510,  1396.913,  1479.978,  1567.982,  1661.219,  1760.000,  1864.655,  1975.533,
  2093.004,  2217.461,  2349.318,  2489.016,  2637.021,  2793.826,  2959.955,  3135.964,  3322.437,  3520.000,  3729.309,  3951.067,
};

gint nharnmonics[G_N_ELEMENTS(midi2frq)];

// Which offset we subtract from the sequence of complex numbers
typedef enum
{
  CENTER_MODE_FIRST = 0,  // No offset compensation
  CENTER_MODE_LAST,       // Convergence point (depends on niter inside the set)
  CENTER_MODE_AVERAGE     // Average of the values
} CenterMode;

// fast sine waves
typedef struct _fastsin
{
  gdouble si0, si1, fc;
} fsin;

typedef struct _complexd {
  gdouble r;
  gdouble i;
} complexd;

typedef struct _AppData
{
  GtkWidget *window;
  guint w, h, y, h2;
  // rendering idle handler
  guint render_id;
  // complex plane
  gdouble crx, crw, crs, cr;
  gdouble ciy, cih, cis, ci;
  // last orbit
  gint niter;

  // track motion?
  gboolean motion;
  // how to map the orbit to the harmonics
  CenterMode center_mode;

  // render buffer
  cairo_surface_t *pix;

  // fractal harmonics settings
  gint ntime, nfreq;
  complexd *v;  // complext orbit path
  complexd *f;  // magnitude and phase

  // gstreamer
  GstElement *pipe;
  GstElement *src;

  // osc for additive synth
  gint oct;   // 0...6
  gdouble *wave;
  // TODO: do per voice
  gint note;  // 0...11
  fsin *fs;

  // osc for harmonics display
  gdouble *waved;
  fsin *fsd;
} AppData;
static AppData app = { 0, };

static guint
do_mandelbrot_traced (gdouble cr, gdouble ci, guint maxn, complexd * v)
{
  // z = z² + c
  gdouble zr = cr, zi = ci, zt;
  guint n = 0;

  while (n < maxn) {
    v[n].r = zr;
    v[n].i = zi;
    zt = cr + zr * zr - zi * zi;
    zi = ci + 2.0 * zi * zr;
    zr = zt;

    if (sqrt (zr * zr + zi * zi) > 2.0) {
      guint m;

      // set the rest to inital position to avoid noise
      for (m = n; m < maxn; m++) {
        v[m].r = zr;
        v[m].i = zi;
      }
      break;
    }
    n++;
  }
  return n;
}


static guint
do_mandelbrot (gdouble cr, gdouble ci, guint maxn)
{
  // z = z² + c
  gdouble zr = cr, zi = ci, zt;
  guint n = 0;

  while (n < maxn) {
    zt = cr + zr * zr - zi * zi;
    zi = ci + 2.0 * zi * zr;
    zr = zt;

    if (sqrt (zr * zr + zi * zi) > 2.0)
      break;
    n++;
  }
  return n;
}

static gboolean
do_mandelbot_image (gpointer user_data)
{
  AppData *self = (AppData *) user_data;
  guint8 *src = cairo_image_surface_get_data (self->pix), *line1, *line2;
  guint stride = cairo_image_surface_get_stride (self->pix);
  guint x, w = self->w, y = self->y, h = self->h, h2 = self->h2;
  guint c;
  gdouble cr = self->crx, crs = self->crs;
  gdouble ci = self->ci;

  // fill one line and do y-mirroring
  line1 = src + (y * stride);
  line2 = src + (((h - 1) - y) * stride);
  for (x = 0; x < w; x++) {
    c = do_mandelbrot (cr, ci, 256);
    if (c == 256) {
      (*line1++) = (*line2++) = 0;
      (*line1++) = (*line2++) = 0;
      (*line1++) = (*line2++) = 0;
    } else {
      (*line1++) = (*line2++) = 64 + (c >> 1);
      (*line1++) = (*line2++) = c >> 1;
      (*line1++) = (*line2++) = c;
    }
    (*line1++) = (*line2++) = 0;
    cr += crs;
  }

  // refresh every 8 lines
  if ((y & 0x7) == 0x7) {
    gtk_widget_queue_draw_area (self->window, 0, y - 7, w, 8);
    gtk_widget_queue_draw_area (self->window, 0, (h - y), w, 8);
  }
  self->y = ++y;
  self->ci += self->cis;
  if (y < h2) {
    return TRUE;
  } else {
    // final refresh
    gtk_widget_queue_draw_area (self->window, 0, y - 7, w, 16);
    return FALSE;
  }
}

static void
process_orbit (AppData *self)
{
  gdouble cr = self->cr, ci = self->ci;
  complexd *v = self->v;
  complexd *f = self->f;
  gdouble tr, ti;
  gdouble mir, mar, mii, mai;
  gint i,j;
  gint niter = self->niter, nfreq = self->nfreq, ntime = self->ntime;

  switch (self->center_mode) {
    case CENTER_MODE_FIRST:
      // do nothing = center on starting point
      break;
    case CENTER_MODE_LAST:
      // center on last value (the converged value in the set)
      cr = v[niter].r;
      ci = v[niter].i;
      break;
    case CENTER_MODE_AVERAGE:
      // center on average
      mir = mar = v[0].r;
      mii = mai = v[0].i;
      for (i = 1; i < niter; i++) {
        if (v[i].r < mir)
          mir = v[i].r;
        else if (v[i].r > mar)
          mar = v[i].r;
        if (v[i].i < mii)
          mii = v[i].i;
        else if (v[i].i > mai)
          mai = v[i].i;
      }
      cr = mir + (mar - mir) / 2.0;
      ci = mii + (mai - mii) / 2.0;
      break;
    default:
      break;
  }
  // convert orbit series to magnitude and phase pairs
  // should we set [0] to some special value (dc offset)?
  for (i = 0; i < niter; i++) {
    tr = v[i].r - cr;
    ti = v[i].i - ci;
    f[i].r = sqrt (tr * tr + ti * ti);
    f[i].i = atan2 (ti, tr);
    /*
       f[i].r = tr;
       f[i].i = ti;
     */
  }
  for (i = niter; i < nfreq; i++) {
    f[i].r = f[i].i = 0.0;
  }

  fsin *fs = self->fsd;
  gdouble *w = self->waved;
  gdouble s, m = 0.0;
  for (i = 0; i < ntime;) {
    s = 0.0;
    for (j = 0; j < nfreq; j++) {
      fs[j].si0 = fs[j].fc * fs[j].si1 - fs[j].si0;
      s += fs[j].si0 * f[j].r;
    }
    w[i++] = s;
    if (fabs (s) > m)
      m = fabs (s);

    s = 0.0;
    for (j = 0; j < nfreq; j++) {
      fs[j].si1 = fs[j].fc * fs[j].si0 - fs[j].si1;
      s += fs[j].si1 * f[j].r;
    }
    w[i++] = s;
    if (fabs (s) > m)
      m = fabs (s);
  }
  for (i = 0; i < ntime; i++) {
    w[i] = w[i] / m;
  }
}

// audio processing helpers

static void
init_fsin (fsin *fs, gdouble cycle_length, gint nfreq)
{
 // phase_inc = 2*PI / cycle_samples
  gdouble angle, base =  (2.0 * M_PI) / cycle_length;
  gint j;

  for (j = 0; j < nfreq; j++) {
    angle = (j + 1) * base;
    fs[j].si0 = sin(-angle);
    fs[j].fc = 2.0 * cos(angle);
  }
}

static void
setup_osc (AppData * self)
{
  gint j;
  fsin *fs = self->fs;

  for (j = 0; j < self->nfreq; j++) {
    fs[j].si0 = fs[j].si1 = fs[j].fc = 0.0;
  }

  if (self->note == -1) {
    return;
  }

  gint note = self->oct * 12 + self->note;
  gdouble frq = midi2frq[note];
  // sampling rate is 44100 (SRATE)
  // if freq = 440, cycle_samples = 44100 / 440
  init_fsin (fs, ((gdouble)SRATE) / frq, MIN(self->nfreq, nharnmonics[note]));
}

static void
setup_harmonics ()
{
  gint i,j;
  gdouble frq, nyq = (SRATE / 2.0);

  // compute how many harmonics we can add without aliasing
  for (i = 0; i < G_N_ELEMENTS(midi2frq); i++) {
    frq = midi2frq[i];
    for (j = 1; (j * frq) < nyq; j++) { }
    nharnmonics[i] = j;
    //printf("%3d, %6lf, %3d\n", i, frq, j);
  }
}

// gfx helpers

static void
cairo_show_right_aligned_text_at (cairo_t * cr, gdouble x, gdouble y,
    const gchar *text)
{
  cairo_text_extents_t ext;

  cairo_text_extents (cr, text, &ext);
  cairo_move_to (cr, x - ext.width, y);
  cairo_show_text (cr, text);
}

static gboolean
on_draw (GtkWidget * widget, cairo_t * cr, gpointer user_data)
{
  AppData *self = (AppData *) user_data;
  const gdouble font_size = 10.0;

  // blit fractal
  cairo_set_source_surface (cr, self->pix, 0.0, 0.0);
  cairo_paint (cr);

  cairo_set_font_size (cr, font_size);

  // settings
  {
    gchar text[100];
    gchar *note_names[] = {
      "c-", "c#", "d-", "d#", "e-", "f-", "f#", "g-", "g#", "a-", "a#", "h-",
    };

    cairo_set_source_rgb (cr, 1.0, 1.0, 1.0);
    sprintf (text, "center: %d", self->center_mode);
    cairo_show_right_aligned_text_at (cr, self->w - 5.0, font_size, text);
    if (self->note == -1) {
      sprintf (text, "note: ---");
    } else {
      sprintf (text, "note: %s%1d", note_names[self->note], self->oct + 1);
    }
    cairo_show_right_aligned_text_at (cr, self->w - 5.0, 2 * font_size, text);
    cairo_stroke (cr);
  }

  // draw overlays once we have data
  if (self->niter == -1) {
    return TRUE;
  }

  guint i;

  process_orbit (self);

  // orbit plot
  {
    complexd *v = self->v;
    gdouble x, rx = self->crx, rw = self->w / self->crw;
    gdouble y, iy = self->ciy, ih = self->h / self->cih;

    cairo_set_source_rgb (cr, 1.0, 1.0, 1.0);
    for (i = 0; i < self->nfreq; i++) {
      x = (v[i].r - rx) * rw;
      y = (v[i].i - iy) * ih;
      //cairo_rectangle (cr, x-0.5,y-0.5,1.0,1.0);
      cairo_arc (cr, x, y, 1.0, 0.0, 2 * M_PI);
      cairo_fill (cr);
    }
    cairo_set_line_width (cr, 0.5);
    cairo_set_source_rgba (cr, 1.0, 1.0, 1.0, 0.5);
    x = (v[0].r - rx) * rw;
    y = (v[0].i - iy) * ih;
    cairo_move_to (cr, x, y);
    for (i = 1; i < self->nfreq; i++) {
      x = (v[i].r - rx) * rw;
      y = (v[i].i - iy) * ih;
      cairo_line_to (cr, x, y);
    }
    cairo_stroke (cr);
  }

  // waveform (top-left, 0..w/3, 0..h/5)
  {
    gdouble hs = self->h / 10.0;
    gdouble xs = self->w / (3.0 * self->ntime);
    gdouble *nt = self->waved;

    cairo_set_source_rgba (cr, 1.0, 1.0, 1.0, 0.15);
    cairo_rectangle (cr, 0.0, 0.0, self->w / 3.0, 2 * hs);
    cairo_fill (cr);

    cairo_set_source_rgb (cr, 1.0, 1.0, 1.0);
    cairo_move_to (cr, 3.0, 2 * hs + font_size);
    cairo_show_text (cr, "waveform");
    cairo_stroke (cr);

    cairo_move_to (cr, 0.0, hs);
    for (i = 0; i < self->ntime; i++) {
      cairo_line_to (cr, i * xs, (nt[i] + 1.0) * hs);
    }
    cairo_stroke (cr);
  }

  // spectrogram (bottom-left, 0..w/3, 4*h/5..h)
  {
    gdouble h = self->h, hs = h / 5.0, v;
    gdouble xs = self->w / (3.0 * self->nfreq);
    complexd *f = self->f;

    cairo_set_source_rgba (cr, 1.0, 1.0, 1.0, 0.15);
    cairo_rectangle (cr, 0.0, h - hs, self->w / 3.0, hs);
    cairo_fill (cr);

    cairo_set_source_rgb (cr, 1.0, 1.0, 1.0);
    cairo_move_to (cr, 3.0, h - (hs + (font_size / 2.0)));
    cairo_show_text (cr, "spectrogram (magnitude: red, phase: green)");
    cairo_stroke (cr);

    // magnitude
    cairo_set_source_rgb (cr, 1.0, 0.3, 0.0);
    cairo_move_to (cr, 0.0, h);
    for (i = 0; i < self->nfreq; i++) {
      v = MIN (f[i].r, 1.0);
      cairo_line_to (cr, i * xs, h - (v * hs));
    }
    cairo_stroke (cr);

    // phase
    cairo_set_source_rgb (cr, 0.3, 1.0, 0.0);
    cairo_move_to (cr, 0.0, h);
    for (i = 0; i < self->nfreq; i++) {
      v = (M_PI + f[i].i) / (M_PI + M_PI);
      cairo_line_to (cr, i * xs, h - (v * hs));
    }
    cairo_stroke (cr);
  }

  return TRUE;
}

static void
on_size_allocate (GtkWidget * widget, GtkAllocation * allocation,
    gpointer user_data)
{
  AppData *self = (AppData *) user_data;

  // stop rendering
  if (self->render_id)
    g_source_remove (self->render_id);

  self->w = allocation->width;
  self->h = allocation->height;
  self->h2 = self->h >> 1;

  if (self->pix)
    cairo_surface_destroy (self->pix);
  self->pix = cairo_image_surface_create (CAIRO_FORMAT_RGB24, self->w, self->h);

  self->crx = -2.0;
  self->crw = 2.5;
  self->crs = self->crw / self->w;
  self->ciy = -1.15;
  self->cih = 2.3;
  self->cis = self->cih / self->h;

  // repaint all
  self->ci = self->ciy;
  self->y = 0;
  self->render_id = g_idle_add (do_mandelbot_image, (gpointer) self);
}

static void
update_note (AppData *self, gint note) {
    if (self->note == note)
      return;

    self->note = note;
    setup_osc (self);
    gtk_widget_queue_draw (self->window);
}

static gboolean
on_interaction_event (GtkWidget * widget, GdkEvent * event, gpointer user_data)
{
  AppData *self = (AppData *) user_data;

  switch (event->type) {
    case GDK_BUTTON_PRESS:
      if (event->button.button == 1) {
        self->motion = TRUE;
      }
      break;
    case GDK_BUTTON_RELEASE:
      if (event->button.button == 1) {
        self->motion = FALSE;
      }
      break;
    case GDK_KEY_PRESS:
      switch (event->key.keyval) {
        case GDK_KEY_1:
          if (self->center_mode > 0) {
            self->center_mode--;
            gtk_widget_queue_draw (self->window);
          }
          break;
        case GDK_KEY_2:
          if (self->center_mode < CENTER_MODE_AVERAGE) {
            self->center_mode++;
            gtk_widget_queue_draw (self->window);
          }
          break;
        case GDK_KEY_3:
          if (self->oct > 0) {
            self->oct--;
            setup_osc (self);
            gtk_widget_queue_draw (self->window);
          }
          break;
        case GDK_KEY_4:
          if (self->oct < 6) {
            self->oct++;
            setup_osc (self);
            gtk_widget_queue_draw (self->window);
          }
          break;
        case GDK_KEY_y:
        case GDK_KEY_z:
          update_note (self, 0);
          break;
        case GDK_KEY_s:
          update_note (self, 1);
          break;
        case GDK_KEY_x:
          update_note (self, 2);
          break;
        case GDK_KEY_d:
          update_note (self, 3);
          break;
        case GDK_KEY_c:
          update_note (self, 4);
          break;
        case GDK_KEY_v:
          update_note (self, 5);
          break;
        case GDK_KEY_g:
          update_note (self, 6);
          break;
        case GDK_KEY_b:
          update_note (self, 7);
          break;
        case GDK_KEY_h:
          update_note (self, 8);
          break;
        case GDK_KEY_n:
          update_note (self, 9);
          break;
        case GDK_KEY_j:
          update_note (self, 10);
          break;
        case GDK_KEY_m:
          update_note (self, 11);
          break;
        default:
          break;
      }
      break;
    case GDK_KEY_RELEASE:
      switch (event->key.keyval) {
        case GDK_KEY_y:
        case GDK_KEY_z:
        case GDK_KEY_s:
        case GDK_KEY_x:
        case GDK_KEY_d:
        case GDK_KEY_c:
        case GDK_KEY_v:
        case GDK_KEY_g:
        case GDK_KEY_b:
        case GDK_KEY_h:
        case GDK_KEY_n:
        case GDK_KEY_j:
        case GDK_KEY_m:
          update_note (self, -1);
          break;
        default:
          break;
      }
      break;
    case GDK_MOTION_NOTIFY:
      if (self->motion) {
        // map event->button.x, event->button.y to cr,ci;
        self->cr = self->crx + event->button.x * self->crs;
        self->ci = self->ciy + event->button.y * self->cis;

        // calculate mandelbrot, but store intermediate values
        self->niter = do_mandelbrot_traced (self->cr, self->ci, self->nfreq,
            self->v);

        gtk_widget_queue_draw (self->window);

      }
      break;
    default:
      break;
  }
  return FALSE;
}

static void
on_need_data (GstAppSrc * appsrc, guint length, gpointer user_data)
{
  AppData *self = (AppData *) user_data;
  gsize size = sizeof (gdouble) * self->ntime;

  // length =  4096, size =  3528
  if (length != -1 && length < size)
    return;

  gint nfreq = self->nfreq, ntime = self->ntime;
  gint i,j;
  fsin *fs = self->fs;
  gdouble s, *w = self->wave;
  complexd *f = self->f;
  for (i = 0; i < ntime;) {
    s = 0.0;
    for (j = 0; j < nfreq; j++) {
      fs[j].si0 = fs[j].fc * fs[j].si1 - fs[j].si0;
      s += fs[j].si0 * f[j].r;
    }
    w[i++] = s;

    s = 0.0;
    for (j = 0; j < nfreq; j++) {
      fs[j].si1 = fs[j].fc * fs[j].si0 - fs[j].si1;
      s += fs[j].si1 * f[j].r;
    }
    w[i++] = s;
  }

  GstBuffer *buf = gst_buffer_new_wrapped_full (GST_MEMORY_FLAG_READONLY,
      (gpointer) (self->wave), size, 0, size, NULL, NULL);

  gst_app_src_push_buffer (GST_APP_SRC (self->src), buf);
}

static void
initialize (AppData * self)
{
  GstElement *conv, *tee, *q1, *q2, *wavenc, *fsink, *asink;
  GstCaps *caps;

  // catch first click, to activate drawing
  self->niter = -1;

  // prepare audio settings
  self->ntime = 512;
  self->nfreq = 100;
  self->v = g_new0 (complexd, self->nfreq);
  self->f = g_new0 (complexd, self->nfreq);

  printf ("Using max %d harmonics, generating audio-chunks of %d samples\n",
      self->nfreq, self->ntime);

  // setup osc
  self->oct = 2;
  self->note = -1;
  self->fs = g_new0 (fsin, self->nfreq);
  self->wave = g_new0 (gdouble, self->ntime);
  setup_harmonics ();

  // spectral display
  self->fsd = g_new0 (fsin, self->nfreq);
  self->waved = g_new0 (gdouble, self->ntime);
  init_fsin (self->fsd, self->ntime, self->nfreq);

  // create window and connect to signals
  self->window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  g_signal_connect (G_OBJECT (self->window), "destroy",
      G_CALLBACK (gtk_main_quit), NULL);
  g_signal_connect (G_OBJECT (self->window), "size-allocate",
      G_CALLBACK (on_size_allocate), (gpointer) self);
  g_signal_connect (G_OBJECT (self->window), "draw",
      G_CALLBACK (on_draw), (gpointer) self);
  g_signal_connect (G_OBJECT (self->window), "event",
      G_CALLBACK (on_interaction_event), (gpointer) self);
  gtk_widget_set_size_request (self->window, 600, 450);
  gtk_window_set_title (GTK_WINDOW (self->window),
      "Mandelbrot synthesizer: press left-button and move, press keys (mind your audio volume)");
  gtk_widget_add_events (self->window,
      GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK |
      GDK_BUTTON1_MOTION_MASK);

  // setup gstreamer pipeline
  self->pipe = gst_pipeline_new (NULL);
  self->src = gst_element_factory_make ("appsrc", NULL);
  caps = gst_caps_new_simple ("audio/x-raw",
      "format", G_TYPE_STRING, GST_AUDIO_NE (F64),
      "rate", G_TYPE_INT, SRATE,
      "layout", G_TYPE_STRING, "interleaved",
      "channels", G_TYPE_INT, 1,
      NULL);
  gst_app_src_set_caps (GST_APP_SRC (self->src), caps);
  gst_app_src_set_stream_type (GST_APP_SRC (self->src),
      GST_APP_STREAM_TYPE_STREAM);
  g_signal_connect (G_OBJECT (self->src), "need-data",
      G_CALLBACK (on_need_data), (gpointer) self);
  conv = gst_element_factory_make ("audioconvert", NULL);
  tee = gst_element_factory_make ("tee", NULL);
  q1 = gst_element_factory_make ("queue", NULL);
  asink = gst_element_factory_make ("autoaudiosink", NULL);
  q2 = gst_element_factory_make ("queue", NULL);
  wavenc = gst_element_factory_make ("wavenc", NULL);
  fsink = gst_element_factory_make ("filesink", NULL);
  g_object_set (G_OBJECT (fsink), "location", "mandelsynth3.wav", NULL);
  g_object_set (G_OBJECT (q1), "max-size-buffers", 1, "max-size-bytes", 0,
      "max-size-time", G_GUINT64_CONSTANT (0), "silent", TRUE, NULL);
  g_object_set (G_OBJECT (q2), "max-size-buffers", 1, "max-size-bytes", 0,
      "max-size-time", G_GUINT64_CONSTANT (0), "silent", TRUE, NULL);

  gst_bin_add_many (GST_BIN (self->pipe), self->src, conv, tee, q1, asink, q2,
      wavenc, fsink, NULL);
  gst_element_link_many (self->src, conv, tee, q1, asink, NULL);
  gst_element_link_many (tee, q2, wavenc, fsink, NULL);
  gst_element_set_state (self->pipe, GST_STATE_PLAYING);
}

static void
finalize (AppData * self)
{
  if (self->pipe) {
    gst_app_src_end_of_stream (GST_APP_SRC (self->src));
    gst_element_set_state (self->pipe, GST_STATE_NULL);
    g_object_unref (self->pipe);
  }

  g_free (self->v);
  g_free (self->f);

  g_free (self->fs);
  g_free (self->wave);

  g_free (self->fsd);
  g_free (self->waved);

  if (self->pix)
    cairo_surface_destroy (self->pix);
}

int
main (int argc, char **argv)
{
  AppData *self = &app;

  gtk_init (&argc, &argv);
  gst_init (&argc, &argv);

  initialize (self);

  // show and run
  gtk_widget_show_all (self->window);
  gtk_main ();

  finalize (self);

  return (0);
}
