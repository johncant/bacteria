#include "simple_cell.h"
#include <armadillo>
#include <iostream>

using namespace Bacteria;
using namespace arma;
using namespace std;

// Helpers declarations
static void cairo_draw_bacteria(cairo_t *cr, double x, double y, double width, double length, double angle, double fgi, double bgi);

// Implementation
SimpleCell::SimpleCell(double x, double y, double angle) : x(x), y(y), angle(angle) {
  initialize_wla();
}


SimpleCell::SimpleCell(double x, double y, double angle, double width, double length) : x(x), y(y), angle(angle), width(width), length(length) {
}


SimpleCell::SimpleCell() {
}


void SimpleCell::initialize(double x, double y, double angle) {
  this->x = x;
  this->y = y;
  this->angle = angle;
  initialize_wla();
}


void SimpleCell::initialize_wla() {
  this->width = 1e-6 + 5e-8*randn();
  this->length = 5e-6 + 2.5e-7 * randn();
}


void SimpleCell::draw_onto_cairo(cairo_t *cr, MicroscopeConfig& mc) {
  cairo_draw_bacteria(cr, x*mc.scale, y*mc.scale, width*mc.scale, length*mc.scale, angle, mc.bg_intensity, mc.fg_intensity);
}

// Helpers implementation
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

bool SimpleCell::collision(SimpleCell& other_cell) {
  // Intersect centrelines and ensure minimum distance
  //
  // X0 = t0 M0 + C0
  // X1 = t1 M1 + C1
  //
  // d^2 = (X - X0)T(X - X0)
  //
  //     = (t0 M0 + C0 - t1 M1 - C1)^T (t0 M0 + C0 - t1 M1 - C1)
  //
  //     = (t0 M0 - t1 M1)^T (t0 M0 - t1 M1) + 2 (t0 M0 - t1 M1)^T (C0 - C1) + (C0 - C1)^T (C0 - C1)
  //
  //     = t0^2 M0^T M0 - 2 t0 t1 M0^T M1 + t1^2 M1^T M1 + 2 (t0 M0 - t1 M1)^T (C0 - C1) + (C0 - C1)^T (C0 - C1)
  //
  // partial d(d^2)/dt0 = 2 t0 M0^T M0 - 2 t1 M0^T M1 + 2 M0^T (C0 - C1) = 0
  //
  // partial d(d^2)/dt1 = 2 t1 M1^T M1 - 2 t0 M0^T M1 - 2 M1^T (C0 - C1) = 0
  //
  // [  M0^T M0      -M0^T M1 ] [ t0 ] = -M0^T (C0 - C1)
  // [ -M0^T M1       M1^T M1 ] [ t1 ] =  M1^T (C0 - C1)
  //
  // BigM T = BigK
  //
  // Find t0, t1, clamp into square between -1 and 1. Done!
  //

  Col<double> m0(2), m1(2), c0(2), c1(2);

  c0[0] = this->x;
  c0[1] = this->y;

  c1[0] = other_cell.x;
  c1[1] = other_cell.y;

  m0[0] = 0.5*(this->length - this->width) * cos(this->angle);
  m0[1] = 0.5*(this->length - this->width) * sin(this->angle);

  m1[0] = 0.5*(other_cell.length - other_cell.width) * cos(other_cell.angle);
  m1[1] = 0.5*(other_cell.length - other_cell.width) * sin(other_cell.angle);

  double m0Tm0, m1Tm1, m0Tm1, m0Tc0Mc1, m1Tc0Mc1, c0Mc1Tc0Mc1;

  m0Tm0 = dot(m0, m0);
  m1Tm1 = dot(m1, m1);
  m0Tm1 = dot(m0, m1);
  m0Tc0Mc1 = dot(m0, c0-c1);
  m1Tc0Mc1 = dot(m1, c0-c1);
  c0Mc1Tc0Mc1 = dot(c0-c1, c0-c1);

  Col<double> big_k(2), t0t1;
  Mat<double> big_m;
  double t0, t1;

  big_m <<  m0Tm0 << -m0Tm1 << endr
        << -m0Tm1 <<  m1Tm1 << endr;

  big_k[0] = -m0Tc0Mc1;
  big_k[1] = m1Tc0Mc1;

  t0t1 = solve(big_m, big_k);

  t0 = t0t1[0];
  t1 = t0t1[1];

  double test_t0_d0 = t0, test_t1_d0 = t1;
  double test_t0_d1 = t0, test_t1_d1 = t1;

  // Move to a side of the square
  if (test_t0_d0 < -1) test_t0_d0 = -1;
  if (test_t0_d0 >= 1) test_t0_d0 = 1;

  // Move to a side of the square
  if (test_t1_d1 < -1) test_t1_d1 = -1;
  if (test_t1_d1 >= 1) test_t1_d1 = 1;

  // Adjust t1 accordingly
  // partial d(d^2)/dt0 = 2 t0 M0^T M0 - 2 t1 M0^T M1 + 2 M0^T (C0 - C1) = 0
  // partial d(d^2)/dt1 = 2 t1 M1^T M1 - 2 t0 M0^T M1 - 2 M1^T (C0 - C1) = 0
  test_t1_d0 = (m1Tc0Mc1 + test_t0_d0*m0Tm1) / m1Tm1;
  test_t0_d1 = (test_t1_d1*m0Tm1 - m0Tc0Mc1) / m0Tm0;

  if (t0 >= -1 && t0 < 1 && t1 >= -1 && t1 < 1) {
    // Do nothing. We found the solution using calculus. return true
    return true; // Shortcut
  } else if (test_t1_d0 >= -1 && test_t1_d0 < 1) {
    // One side fixed, one side found using calculus
    t0 = test_t0_d0;
    t1 = test_t1_d0;
  } else if (test_t0_d1 >= -1 && test_t0_d1) {
    // Evaluate distance
    t0 = test_t0_d1;
    t1 = test_t1_d1;
  } else {
    // Use the corner
    if (test_t1_d0 < -1) test_t1_d0 = -1;
    if (test_t1_d0 >= 1) test_t1_d0 = 1;
    t0 = test_t0_d0;
    t1 = test_t1_d0;
  }

  // d = t0^2 M0^T M0 - 2 t0 t1 M0^T M1 +t1^2 M1^T M1 + 2 (t0 M0 - t1 M1)^T (C0 - C1) + (C0 - C1)^T (C0 - C1)

//    double d = t0*t0*m0Tm0 - 2*t0*t1*m0Tm1 + t1*t1*m1Tm1 + 2*t0*m0Tc0Mc1 - 2*t1*m1Tc0Mc1 + c0Mc1Tc0Mc1;
  double d = dot(c0+m0*t0 - c1 - m1*t1, c0+m0*t0 - c1 - m1*t1);
  if (d > 0) {
    d = sqrt(d);
  } else {
    d = 0;
  }

  if (d < 0.5*(this->width + other_cell.width)) {
    return true;
  } else {
    return false;
  }

}
