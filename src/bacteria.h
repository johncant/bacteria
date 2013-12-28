#ifndef __BACTERIA_H__
#define __BACTERIA_H__

#include <armadillo>
#include <cairo.h>
#include "microscope_config.h"

namespace Bacteria {

//  class SimpleCell;
//  template <class CellT> class Colony;
//  template <class CellT> class Generator;
//  class MicroscopeConfig;

  void cairo_draw_background(cairo_t *cr, MicroscopeConfig& mc);

}

#endif
