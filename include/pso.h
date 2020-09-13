#pragma once

#include <stdbool.h>

typedef struct dimension {
    double min;
    double max;
} dimension_t;

typedef struct particle particle_t;
struct particle {
    double *x;                  // an array of nr_dimensions elements

    dimension_t *dimensions;    // an array of nr_dimensions elements
    int nr_dimensions;

    double fitness;
    double *velocity;           // an array of nr_dimensions elements
    particle_t *previous_best;
};

bool particle_create(particle_t **particle, int nr_dimensions);
void particle_destroy(particle_t *particle);
