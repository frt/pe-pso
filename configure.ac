#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_CONFIG_MACRO_DIRS([m4])
AC_INIT([pe-pso], [0.0.0], [fr.teodoro@gmail.com])
AC_CONFIG_SRCDIR([include/pso.h])
AC_LANG(C)
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT

# Checks for programs.
AC_PROG_CC
AC_CHECK_PROGS([MPICC], [mpicc], [:])
if test "$MPICC" = :; then
    AC_MSG_ERROR([This package needs mpicc.])
fi

# Checks for libraries.
AC_CHECK_LIB([m], [cos], , AC_MSG_ERROR([unable to find the cos() function]))
AC_CHECK_LIB([mtwist], [mt_seed], , AC_MSG_ERROR([mtwist library not found]))

AC_CHECK_LIB([parallel_evolution], [parallel_evolution_run], [],
             AC_MSG_ERROR([parallel_evolution library not found]), [-L/usr/lib/openmpi -lmpi])
#AC_CHECK_LIB([topology_parser], [topology_parser_parse])

# Checks for header files.
AC_CHECK_HEADER([parallel_evolution.h],,AC_MSG_ERROR([parallel_evolution.h not found]))
AC_CHECK_HEADER([mtwist.h], [], [AC_MSG_ERROR(mtwist.h missing (get it at https://github.com/frt/mtwist))])
AC_CHECK_HEADER([randistrs.h], [], [AC_MSG_ERROR(randistrs.h missing (get it at https://github.com/frt/mtwist))])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/topology_parser/Makefile])
AC_OUTPUT
