#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_CONFIG_MACRO_DIRS([m4])

AC_PREREQ([2.69])
AC_INIT([bacteria], [0.0.1], [a.johncant@gmail.com])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/generator.cc])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Check for Armadillo
# Thankyou GNU radio!
GR_ARMADILLO(,,AC_MSG_ERROR([Armadillo not found]))

# Gtk
PKG_CHECK_MODULES([gtk], [gtk+-3.0],, [
                      PKG_CHECK_MODULES([gtk], [gtk+-2.0],,AC_MSG_ERROR([gtk+-2.0 or gtk+-3.0 not found]))
                      ])

PKG_CHECK_MODULES([cairo], [cairo])

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
