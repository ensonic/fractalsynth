/*
 * mandelbrot synthesizer
 *
 * Copyright (C) 2011 Stefan Kost <ensonic@users.sf.net>
 *
 * gcc -Wall -g  mandelsynth.c -o mandelsynth `pkg-config gtk+-2.0 cairo gstreamer-0.10 gstreamer-app-0.10 gstreamer-fft-0.10 --cflags --libs` -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>
#include <gtk/gtk.h>
#include <gst/gst.h>
#include <gst/app/gstappsrc.h>
#include <gst/fft/gstfftf64.h>

typedef struct _AppData
{
  GtkWidget *window;
  guint w, h, y, h2;
  // rendering idle handler
  guint render_id;
  // complex plane
  gdouble crx, crw, crs;
  gdouble ciy, cih, cis, ci;

  // track motion?
  gboolean motion;

  // render buffer
  cairo_surface_t *pix;

  // fft
  gint ntime, nfreq;
  GstFFTF64 *ifft;
  GstFFTF64Complex *v;
  GstFFTF64Complex *f;
  gdouble *t, *nt;

  // gstreamer
  GstElement *pipe;
  GstElement *src;

} AppData;
static AppData app = { 0, };

static guint
do_mandelbrot_traced (gdouble cr, gdouble ci, guint maxn, GstFFTF64Complex * v)
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

      // set the rest to 0.0 to avoid noise
      for (m = n; m < maxn; m++) {
        v[m].r = v[m].i = 0.0;
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
render_fractal_image (gpointer user_data)
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

static gboolean
on_expose (GtkWidget * widget, GdkEventExpose * event, gpointer user_data)
{
  AppData *self = (AppData *) user_data;
  cairo_t *cr;

  if (!gtk_widget_get_realized (self->window))
    return TRUE;

  cr = gdk_cairo_create (event->window);

  // blit fractal
  gdk_cairo_region (cr, event->region);
  cairo_clip (cr);

  cairo_set_source_surface (cr, self->pix, 0.0, 0.0);
  cairo_paint (cr);

  // draw overlays
  if (self->motion) {
    guint i;

    // orbit plot
    {
      GstFFTF64Complex *v = self->v;
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
      gdouble *nt = self->nt;

      cairo_set_source_rgba (cr, 1.0, 1.0, 1.0, 0.15);
      cairo_rectangle (cr, 0.0, 0.0, self->w / 3.0, 2 * hs);
      cairo_fill (cr);

      cairo_set_source_rgb (cr, 1.0, 1.0, 1.0);
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
      GstFFTF64Complex *f = self->f;

      cairo_set_source_rgba (cr, 1.0, 1.0, 1.0, 0.15);
      cairo_rectangle (cr, 0.0, h - hs, self->w / 3.0, hs);
      cairo_fill (cr);

      // magnitude
      cairo_set_source_rgb (cr, 1.0, 0.5, 0.0);
      cairo_move_to (cr, 0.0, h);
      for (i = 0; i < self->nfreq; i++) {
        v = MIN (f[i].r, 1.0);
        cairo_line_to (cr, i * xs, h - (v * hs));
      }
      cairo_stroke (cr);

      // phase
      cairo_set_source_rgb (cr, 0.5, 1.0, 0.0);
      cairo_move_to (cr, 0.0, h);
      for (i = 0; i < self->nfreq; i++) {
        v = (M_PI + f[i].i) / (M_PI + M_PI);
        cairo_line_to (cr, i * xs, h - (v * hs));
      }
      cairo_stroke (cr);
    }
  }

  cairo_destroy (cr);
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
  self->render_id = g_idle_add (render_fractal_image, (gpointer) self);
}

static gboolean
on_pointer_event (GtkWidget * widget, GdkEvent * event, gpointer user_data)
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
    case GDK_MOTION_NOTIFY:
      if (self->motion) {
        // map event->button.x, event->button.y to cr,ci;
        gdouble cr = self->crx + event->button.x * self->crs;
        gdouble ci = self->ciy + event->button.y * self->cis;
        GstFFTF64Complex *v = self->v;
        GstFFTF64Complex *f = self->f;
        gdouble *t = self->t;
        gdouble *nt = self->nt;
        gdouble tr, ti;
        gdouble m = 0.0;
        gdouble mir, mar, mii, mai;
        gint i;

        // calculate mandelbrot, but store intermediate values
        do_mandelbrot_traced (cr, ci, self->nfreq, v);
        gtk_widget_queue_draw (self->window);
        // determine center (this is not necessarily close to cr,ci
        mir = mar = v[0].r;
        mii = mai = v[0].i;
        for (i = 1; i < self->nfreq; i++) {
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
        // convert orbit series to magnitude and phase pairs
        // should we set the 1st bin to some special value (dc offset)
        for (i = 0; i < self->nfreq; i++) {
          tr = v[i].r - cr;
          ti = v[i].i - ci;
          f[i].r = sqrt (tr * tr + ti * ti);
          f[i].i = atan2 (ti, tr);
          /*
             f[i].r = tr;
             f[i].i = ti;
           */
        }
        // do iFFT of the data
        gst_fft_f64_inverse_fft (self->ifft, f, t);
        // normalize
        for (i = 0; i < self->ntime; i++) {
          if (fabs (t[i]) > m)
            m = fabs (t[i]);
        }
        for (i = 0; i < self->ntime; i++) {
          nt[i] = t[i] / m;
        }
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
  GstBuffer *buf;
  guint size;

  size = sizeof (gdouble) * self->ntime;

  if (length != -1 && length < size)
    return;

  buf = gst_buffer_new ();
  GST_BUFFER_DATA (buf) = (guint8 *) (self->nt);
  GST_BUFFER_SIZE (buf) = size;

  gst_app_src_push_buffer (GST_APP_SRC (self->src), buf);
}

static void
initialize (AppData * self)
{
  GstElement *conv, *tee, *q1, *q2, *wavenc, *fsink, *asink;
  GstCaps *caps;

  // prepare fft
  self->ntime = gst_fft_next_fast_length (250);
  self->nfreq = self->ntime / 2 + 1;
  self->ifft = gst_fft_f64_new (self->ntime, TRUE);
  self->v = g_new0 (GstFFTF64Complex, self->nfreq);
  self->f = g_new0 (GstFFTF64Complex, self->nfreq);
  self->t = g_new0 (gdouble, self->ntime);
  self->nt = g_new0 (gdouble, self->ntime);

  printf ("Using a FFT with %d bands, generating audio-chunks of %d samples\n",
      self->nfreq, self->ntime);

  // create window and connect to signals
  self->window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  g_signal_connect (G_OBJECT (self->window), "destroy",
      G_CALLBACK (gtk_main_quit), NULL);
  g_signal_connect (G_OBJECT (self->window), "size-allocate",
      G_CALLBACK (on_size_allocate), (gpointer) self);
  g_signal_connect (G_OBJECT (self->window), "expose-event",
      G_CALLBACK (on_expose), (gpointer) self);
  g_signal_connect (G_OBJECT (self->window), "event",
      G_CALLBACK (on_pointer_event), (gpointer) self);
  gtk_widget_set_size_request (self->window, 600, 450);
  gtk_window_set_title (GTK_WINDOW (self->window),
      "Mandelbrot synthesizer: press left-button over black set and move (mind your audio volume)");
  gtk_widget_add_events (self->window,
      GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK |
      GDK_BUTTON1_MOTION_MASK);

  // setup gstreamer pipeline
  self->pipe = gst_pipeline_new (NULL);
  self->src = gst_element_factory_make ("appsrc", NULL);
  caps = gst_caps_new_simple ("audio/x-raw-float",
      "rate", G_TYPE_INT, 44100,
      "channels", G_TYPE_INT, 1,
      "endianness", G_TYPE_INT, BYTE_ORDER, "width", G_TYPE_INT, 64, NULL);
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
  g_object_set (G_OBJECT (fsink), "location", "mandelsynth.wav", NULL);

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

  if (self->ifft)
    gst_fft_f64_free (self->ifft);
  g_free (self->v);
  g_free (self->f);
  g_free (self->t);
  g_free (self->nt);

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
