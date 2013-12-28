#include "bacteria.h"

// Implementation
void Bacteria::cairo_draw_background(cairo_t *cr, MicroscopeConfig& mc) {

  cairo_pattern_t *bg = cairo_pattern_create_rgb(mc.bg_intensity, mc.bg_intensity, mc.bg_intensity);
  cairo_rectangle(cr, 0, 0, mc.view_width, mc.view_height);
  cairo_set_source(cr, bg);

  cairo_fill(cr);
}

