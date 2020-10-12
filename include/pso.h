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

struct particle {
    double *x;                  // an array of nr_dimensions elements
    double fitness;
    double *velocity;           // an array of nr_dimensions elements

    double *previous_best_x;      // previous best position
    double previous_best_fitness;

    // particles informed by this one
    particle_t **neighbours;    // an array of pointers to particles

    double *neighbourhood_best_x;
    double neighbourhood_best_fitness;
};

typedef struct pso_config {
    limits_t *search_space_limits;
    int nr_particles;
    int nr_neighbours;
} pso_config_t;

typedef struct swarm {
    particle_t **particles;
    int nr_particles;
    limits_t *search_space;
    double best_fitness;
    int best_fitness_iteration;
} swarm_t;

/**
 * Creates an empty new particle.
 *
 * \return The pointer to the particle created or NULL if creation fails.
 */
particle_t *particle_create(pso_config_t *pso_config);
void particle_destroy(particle_t *particle);

/**
 * Creates the initial swarm.
 *
 * \return The pointer to the swarm created or NULL if failed.
 */
swarm_t *swarm_create(pso_config_t *pso_config, double (*fitness_func)(double *x));
void swarm_destroy(swarm_t *swarm);

/**
 * Run iterations.
 *
 * \param swarm The swarm
 * \param pso_config PSO configurations
 * \param nr_iterations The maximum number of iterations
 *
 * \return The number of iterations that where run.
 */
int iterarions(swarm_t *swarm, pso_config_t *pso_config, int nr_iterations);
