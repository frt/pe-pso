#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([pe-pso], [0.0.0], [fr.teodoro@gmail.com])
AC_CONFIG_SRCDIR([include/pso.h])
AC_LANG(C)
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Checks for programs.
AC_PROG_CC
AC_CHECK_PROGS([MPICC], [mpicc], [:])
if test "$MPICC" = :; then
    AC_MSG_ERROR([This package needs mpicc.])
fi

# Checks for libraries.
AC_CHECK_LIB([m], [cos], , AC_MSG_ERROR([unable to find the cos() function]))
AC_CHECK_LIB([parallel_evolution], [parallel_evolution_run], [],
             AC_MSG_ERROR([parallel_evolution library not found]), [-L/usr/lib/openmpi -lmpi])
#AC_CHECK_LIB([topology_parser], [topology_parser_parse])

# Checks for header files.
AC_CHECK_HEADER([parallel_evolution.h],,AC_MSG_ERROR([parallel_evolution.h not found]))

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
