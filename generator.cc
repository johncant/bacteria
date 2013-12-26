#include "generator.h"
#include <iostream>
#include <gtk/gtk.h>

using namespace Bacteria;
using namespace arma;
using namespace std;

static gboolean on_cairo_draw_event(GtkWidget *widget, cairo_t *cr, gpointer user_data);
static void cairo_draw_bacteria(cairo_t *cr, double x, double y, double width, double length, double angle, double fgi, double bgi);
static void cairo_draw_background(cairo_t *cr, double width, double length, double bgi);

Colony<SimpleCell> *test_cells;

//SimpleCell
SimpleCell::SimpleCell(double x, double y, double angle) : x(x), y(y), angle(angle) {
  this->width = 1e-6 + 1e-7*randn();
  this->length = 5e-6 + 2e-6 * randn();
}

void SimpleCell::draw_onto_cairo(cairo_t *cr, MicroscopeConfig& mc) {
  cairo_draw_bacteria(cr, x*mc.scale, y*mc.scale, width*mc.scale, length*mc.scale, angle, mc.fg_intensity, mc.bg_intensity);
}

//Colony

template <class CellT> Colony<CellT> Colony<CellT>::generate_gaussian(double mu=0e0, double sigma=4e-6, int count=20) {
  cout << "Not yet implemented" << endl;
}

template <class CellT> void Colony<CellT>::populate_gaussian(double mu=0e0, double sigma=4e-6, int count=20) {
  for(int i=0; i<count; i++) {
    double cx = mu+sigma*randn(),
           cy = mu+sigma*randn();
    this->push_back(CellT(cx, cy, randu()*2*M_PI));
  }
}

template <class CellT> void Colony<CellT>::draw_onto_cairo(cairo_t* cr, MicroscopeConfig& mc) {
  for(typename vector<CellT>::iterator i = this->begin(); i != this->end(); ++i) {
    i->draw_onto_cairo(cr, mc);
  }
}

//Generator
//
//
//#template <class CellT> Mat<double> Generator<CellT>::generate_image(int sx=400, int sy=400) {
//#  Mat<double> img(sx, sy, fill::zeros);
//#  Colony<CellT> cells;
//#  cells.populate_gaussian(0e0, 4e-6, 100);
//#//  cells.cairo_draw_onto();
//#  MicroscopeConfig mc;
//#  mc.scale = 4e-5;
//#  mc.bg_intensity = 0.6;
//#  mc.fg_intensity = 0.4;
//#  s
//#  return Mat<double>(img);
//#}

template <class CellT> Colony<CellT> Generator<CellT>::generate(double mu, double sigma, int n) {
  Colony<CellT> cells;
  cells.populate_gaussian(mu, sigma, n);
  return Colony<CellT>(cells);
}

template <class CellT> void Generator<CellT>::draw_onto_cairo(cairo_t* cr, Colony<CellT> &cells) {
  MicroscopeConfig mc;
  mc.scale = 1e7;
  mc.bg_intensity = 0.6;
  mc.fg_intensity = 0.4;

  cairo_draw_background(cr, 400, 400, 0.4);
  cells.draw_onto_cairo(cr, mc);
}

// Callbacks

static gboolean on_cairo_draw_event(GtkWidget *widget, cairo_t *cr, gpointer user_data) {

  cairo_draw_background(cr, 400, 400, 0.4);
  Generator<SimpleCell>::draw_onto_cairo(cr, *test_cells);
//  cairo_draw_bacteria(cr, 0, 0, 30, 120, M_PI/6, 0.4, 0.6);
//  cairo_draw_bacteria(cr, -100, 0, 30, 120, M_PI/2, 0.4, 0.6);
//  cairo_draw_bacteria(cr, 0, -100, 30, 120, 0, 0.4, 0.6);

//  cairo_draw_background(cr, 400, 400, 0.4);
//  cairo_draw_bacteria(cr, 0, 0, 30, 120, M_PI/6, 0.4, 0.6);
//  cairo_draw_bacteria(cr, -100, 0, 30, 120, M_PI/2, 0.4, 0.6);
//  cairo_draw_bacteria(cr, 0, -100, 30, 120, 0, 0.4, 0.6);
  return FALSE;
}

// Helpers

static void cairo_draw_background(cairo_t *cr, double width, double length, double bgi) {
  cairo_identity_matrix(cr);
  cairo_translate(cr, 200, 200);

  cairo_pattern_t *bg = cairo_pattern_create_rgb(bgi, bgi, bgi);
  cairo_rectangle(cr, -width/2, -length/2, width, length);
  cairo_set_source(cr, bg);
  cairo_fill(cr);
}

static void cairo_draw_bacteria(cairo_t *cr, double x, double y, double width, double length, double angle, double bgi, double fgi) {

  // x, y, w, l, angle
  //
  cout << "Drawing at" << x << ", " << y << endl;

  // Save the transformation matrix
  cairo_matrix_t original_matrix;
  cairo_get_matrix(cr, &original_matrix);

  double radius = width*0.5;
  double recthlength = (length-width)*0.5;

  cairo_pattern_t *radpattern = cairo_pattern_create_radial(-recthlength, 0, 0, -recthlength, 0, radius);
  cairo_pattern_t *radpattern1 = cairo_pattern_create_radial(recthlength, 0, 0, recthlength, 0, radius);
  cairo_pattern_t *rectpattern = cairo_pattern_create_linear(-radius, -radius, -radius, radius);

  cairo_translate(cr, x, -y);
  cairo_rotate(cr, -angle);

  double j;

  for ( j = -1; j <= 1; j += 0.1 ) {
    double cj = fgi + (bgi-fgi)*cos(j*0.5*M_PI);
    cairo_pattern_add_color_stop_rgb(rectpattern, 0.5+0.5*sin(j*0.5*M_PI), cj, cj, cj);
  }

  for ( j = 0; j <= 1; j += 0.1 ) {
    double f = (bgi-fgi);
    double cj = fgi+f*cos(j*0.5*M_PI);
    double along = sin(j*0.5*M_PI);
    cairo_pattern_add_color_stop_rgb(radpattern, along, cj, cj, cj);
    cairo_pattern_add_color_stop_rgb(radpattern1, along, cj, cj, cj);
  }


  // Rectangle
  cairo_rectangle(cr, -recthlength, -radius, 2*recthlength, width);
  cairo_set_source(cr, rectpattern);
  cairo_fill(cr);
  // Arcs
  cairo_arc(cr, -recthlength, 0, radius, M_PI/2, 3*M_PI/2);
  cairo_set_source(cr, radpattern);
  cairo_fill(cr);
  // Arcs
  cairo_arc_negative(cr, recthlength, 0, radius, M_PI/2, 3*M_PI/2);
  cairo_set_source(cr, radpattern1);
  cairo_fill(cr);

  // Origin
//  cairo_pattern_t *fg = cairo_pattern_create_rgb(1, 0, 0);
//  cairo_rectangle(cr, -5, -5, 10, 10);
//  cairo_set_source(cr, fg);
//  cairo_fill(cr);

  cairo_set_matrix(cr, &original_matrix);

  cairo_pattern_destroy(rectpattern);
  cairo_pattern_destroy(radpattern);
  cairo_pattern_destroy(radpattern1);
}

void test_with_gtk(int argc, char** argv) {
  GtkWidget *window;
  GtkWidget *darea; 

  gtk_init(&argc, &argv);
  window = gtk_window_new(GTK_WINDOW_TOPLEVEL);

  darea = gtk_drawing_area_new();
  gtk_container_add(GTK_CONTAINER (window), darea);

  g_signal_connect(G_OBJECT(darea), "draw", 
            G_CALLBACK(on_cairo_draw_event), NULL);  
  g_signal_connect(G_OBJECT(window), "destroy",
            G_CALLBACK(gtk_main_quit), NULL);

  gtk_widget_set_app_paintable(window, TRUE);

  gtk_widget_show_all(window);

  gtk_window_resize(GTK_WINDOW(window), 400, 400);

  gtk_main();

}


// Test

#include <cstdlib>
#include <time.h>

int main(int argc, char** argv) {
  srand(time(NULL));
  Generator<SimpleCell> gen;
  Colony<SimpleCell> cells = gen.generate(0e0, 5e-6, 20);
  test_cells = &cells;

  test_with_gtk(argc, argv);

  return 0;
}
