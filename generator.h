#include <armadillo>
#include <vector>
#include <cairo.h>


namespace Bacteria {

  class SimpleCell;
  template <class CellT> class Colony;
  template <class CellT> class Generator;
  class MicroscopeConfig;

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
    private:
      void initialize_wla();
  };

  template <class CellT> class Colony : public std::vector<CellT> {
    public:
      static Colony<CellT> generate_gaussian(double mu, double sigma, int count);
      void populate_gaussian(double mu, double sigma, int count);
      void draw_onto_cairo(cairo_t *ct, MicroscopeConfig& mc);
      bool collision(CellT &new_cell);
  };

  template <class CellT> class Generator {
    public:
      static arma::Mat<double> generate_image(int sx, int sy);
      static Colony<CellT> generate(double mu, double sigma, int count);
      static void draw_onto_cairo(cairo_t *ct, Colony<CellT> &cells);
  };

  class MicroscopeConfig {
    public:
      double scale;
      double fg_intensity, bg_intensity;
  };
}


