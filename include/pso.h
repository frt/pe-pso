#pragma once

typedef struct particle particle_t;

typedef struct dimension {
    double min;
    double max;
} dimension_t;

typedef struct limits {
    dimension_t *dimensions;    // an array of nr_dimensions elements
    int nr_dimensions;
} limits_t;

typedef struct neighbourhood {
    particle_t **neighbours;    // an array of pointers to particles
    int nr_neighbours;
} neighbourhood_t;

struct particle {
    double *x;                  // an array of nr_dimensions elements
    double fitness;
    double *velocity;           // an array of nr_dimensions elements
    particle_t *previous_best;

    neighbourhood_t neighbours;
};

typedef struct swarm {
    particle_t *particles;
    int nr_particles;
    limits_t search_space;
} swarm_t;

/**
 * Creates an empty new particle.
 *
 * \return The pointer to the particle created or NULL if creation fails.
 */
particle_t *particle_create(limits_t *dimensions_limits);

void particle_destroy(particle_t *particle);
