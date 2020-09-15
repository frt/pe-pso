#pragma once

#include <stdbool.h>

typedef struct dimension {
    double min;
    double max;
} dimension_t;

typedef struct limits {
    dimension_t *dimensions;    // an array of nr_dimensions elements
    int nr_dimensions;
} limits_t;

typedef struct particle particle_t;
struct particle {
    double *x;                  // an array of nr_dimensions elements

    double fitness;
    double *velocity;           // an array of nr_dimensions elements
    particle_t *previous_best;
};

/**
 * Creates an empty new particle.
 */
particle_t *particle_create(limits_t *dimensions_limits);

void particle_destroy(particle_t *particle);
