AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS=gtk_example

gtk_example_SOURCES=src/simple_cell.cc\
										src/bacteria.cc\
										src/colony.cc\
										src/generator.cc\
										src/gtk_example.cc

AM_CXXFLAGS = $(cairo_CFLAGS) $(gtk_CFLAGS) $(ARMADILLO_CFLAGS)
LIBS = $(cairo_LIBS) $(gtk_LIBS) $(ARMADILLO_LIBS)
