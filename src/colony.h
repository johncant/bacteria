#ifndef __BACTERIA__COLONY_H__
#define __BACTERIA__COLONY_H__

#include "microscope_config.h"
#include "colony.h"
#include <iostream>
#include <vector>
#include <cairo.h>
#include <armadillo>


using namespace Bacteria;
using namespace arma;
using namespace std;

namespace Bacteria {

  template <class CellT> class Colony : public std::vector<CellT> {
    public:
      static Colony<CellT> generate_gaussian(double mu, double sigma, int count);
      void populate_gaussian(double mu, double sigma, int count);
      void draw_onto_cairo(cairo_t *ct, MicroscopeConfig& mc);
      bool collision(CellT &new_cell);
  };

  template <class CellT> Colony<CellT> Colony<CellT>::generate_gaussian(double mu=0e0, double sigma=4e-6, int count=20) {
    std::cout << "Not yet implemented" << std::endl;
  }

  template <class CellT> void Colony<CellT>::populate_gaussian(double mu=0e0, double sigma=4e-6, int count=20) {
    for(int i=0; i<count; i++) {

      CellT new_cell;

      do {
        std::cout << "placing bacteria " << i << std::endl;
        double cx = mu+sigma*arma::randn(),
               cy = mu+sigma*arma::randn();
        new_cell.initialize(cx, cy, arma::randu()*2*M_PI);
      } while (collision(new_cell));

      this->push_back(new_cell);
    }
  }

  template <class CellT> bool Colony<CellT>::collision(CellT &new_cell) {

    for(typename vector<CellT>::iterator i = this->begin(); i != this->end(); ++i) {
      if (new_cell.collision(*i)) {
        return true;
      }
    }

    return false;
  }

  template <class CellT> void Colony<CellT>::draw_onto_cairo(cairo_t* cr, MicroscopeConfig& mc) {
    for(typename vector<CellT>::iterator i = this->begin(); i != this->end(); ++i) {
      i->draw_onto_cairo(cr, mc);
    }
  }

}

#endif
