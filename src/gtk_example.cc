#include <iostream>
#include <time.h>
#include <gtk/gtk.h>

#include "bacteria.h"
#include "generator.h"
#include "colony.h"
#include "simple_cell.h"

using namespace std;
using namespace Bacteria;


// Callback/helper declarations
static gboolean on_cairo_draw_event(GtkWidget *widget, cairo_t *cr, gpointer user_data);
void test_with_gtk(int argc, char** argv);


// Globals for test purposes :-O
Colony<SimpleCell> *test_cells;
MicroscopeConfig *test_mc;

// Callbacks

static gboolean on_cairo_draw_event(GtkWidget *widget, cairo_t *cr, gpointer user_data) {

  cairo_draw_background(cr, *test_mc);
  cairo_translate(cr, test_mc->view_width/2, test_mc->view_height/2);
  test_cells->draw_onto_cairo(cr, *test_mc);

  return FALSE;
}

// Helpers

void test_with_gtk(int argc, char** argv) {
  GtkWidget *window;
  GtkWidget *darea; 

  gtk_init(&argc, &argv);
  window = gtk_window_new(GTK_WINDOW_TOPLEVEL);

  darea = gtk_drawing_area_new();
  gtk_container_add(GTK_CONTAINER (window), darea);

#if GTK_MAJOR_VERSION == 2
  g_signal_connect(G_OBJECT(darea), "expose-event", 
            G_CALLBACK(on_cairo_draw_event), NULL);  
#elif GTK_MAJOR_VERSION == 3
  g_signal_connect(G_OBJECT(darea), "draw", 
            G_CALLBACK(on_cairo_draw_event), NULL);  
#endif

  g_signal_connect(G_OBJECT(window), "destroy",
            G_CALLBACK(gtk_main_quit), NULL);

  gtk_widget_set_app_paintable(window, TRUE);

  gtk_widget_show_all(window);

  gtk_window_resize(GTK_WINDOW(window), 400, 400);

  gtk_main();

}

// ..
int main(int argc, char** argv) {

  // Random seed
  int seed = time(NULL);
  cout << "Random seed: " << seed << endl;
  srand(seed);

  // Generate the bacteria
  Generator<SimpleCell> gen;
  Colony<SimpleCell> cells = gen.generate(0e0, 5e-6, 50);

  // Configure microscope
  MicroscopeConfig mc;
  mc.scale = 1e7;
  mc.bg_intensity = 0.4;
  mc.fg_intensity = 0.6;
  mc.view_width = 400;
  mc.view_height = 400;

  test_mc = &mc;
  test_cells = &cells;

  test_with_gtk(argc, argv);

  return 0;
}


