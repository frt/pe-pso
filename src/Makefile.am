CC=$(MPICC)
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(prefix)/include
AM_CFLAGS = -std=gnu89	# need this for compatibility with libmtwist

lib_LTLIBRARIES = libpso.la
libpso_la_SOURCES = pso.c
nobase_include_HEADERS = ../include/pso.h

bin_PROGRAMS = pe-pso
pe_pso_SOURCES = main.c
pe_pso_LDADD = libpso.la \
			   -L$(libdir) \
			   -lparallel_evolution \
			   -lmtwist \
			   -lconfig \
			   -lm
