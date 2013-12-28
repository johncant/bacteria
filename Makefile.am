AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS=generator
generator_SOURCES=src/generator.cc

AM_CXXFLAGS = $(cairo_CFLAGS) $(gtk_CFLAGS) $(ARMADILLO_CFLAGS)
LIBS = $(cairo_LIBS) $(gtk_LIBS) $(ARMADILLO_LIBS)
