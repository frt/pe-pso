# Particle Swarm Optimization (PSO) using the Parallel Evolution library

This is an implementation of the
[Standard PSO 2011 (SPSO-2011)](http://clerc.maurice.free.fr/pso/SPSO_descriptions.pdf) for use in an
ecosystem of algorithms running in parallel.  The
[parallel_evolution-lib](http://github.com/frt/parallel_evolution-lib) provides
utility functions for running the evolutionary algorithms in parallel and
exchanging candidate solutions between them, besides that it defines how and
when the algorithms communicate with each other and an unified way of reading
configurations, log levels and reporting results.

## Build and Install

### Dependencies

 * [parallel_evolution-lib](http://github.com/frt/parallel_evolution-lib)

### Generate the './configure' script

```bash
aclocal
autoconf
automake --add-missing
```

### Configure and install

```bash
./configure
make
sudo make install
```
