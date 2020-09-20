CC=$(MPICC)
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(prefix)/include
bin_PROGRAMS = pe-pso
pe_pso_SOURCES = pso.c main.c
pe_pso_LDADD = -L$(libdir) \
			   -lparallel_evolution \
			   -ltopology_parser \
			   -lmtwist \
			   -lm
