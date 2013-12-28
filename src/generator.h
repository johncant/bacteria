#ifndef __BACTERIA__GENERATOR_H__
#define __BACTERIA__GENERATOR_H__

#include <armadillo>
#include <cairo.h>
#include "colony.h"

namespace Bacteria {

  template <class CellT> class Generator {
    public:
      static arma::Mat<double> generate_image(int sx, int sy);
      static Colony<CellT> generate(double mu, double sigma, int count);
  };

  template <class CellT> Colony<CellT> Generator<CellT>::generate(double mu, double sigma, int n) {
    Colony<CellT> cells;
    cells.populate_gaussian(mu, sigma, n);
    return Colony<CellT>(cells);
  }

}

#endif
