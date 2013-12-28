#ifndef __BACTERIA__SIMPLE_CELL_H__
#define __BACTERIA__SIMPLE_CELL_H__

#include "microscope_config.h"
#include <cairo.h>

namespace Bacteria {

  class SimpleCell {
    public:
      double width;
      double length;
      double x;
      double y;
      double angle;
      SimpleCell(double x, double y, double angle, double width, double length);
      SimpleCell(double x, double y, double angle);
      SimpleCell();
      void draw_onto_cairo(cairo_t *ct, MicroscopeConfig& mc);
      void initialize(double x, double y, double angle);
      bool collision(SimpleCell &other_cell);
    private:
      void initialize_wla();
  };

}

#endif
