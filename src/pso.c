#include "pso.h"

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include <randistrs.h>

particle_t *particle_create(pso_config_t *pso_config)
{
    particle_t *new;
    int nr_dimensions = pso_config->search_space_limits->nr_dimensions;
    int nr_neighbours = pso_config->nr_neighbours;

    new = (particle_t *)calloc(1, sizeof(particle_t));
    if (new == NULL)
        return NULL;

    new->x = (double *)malloc(nr_dimensions * sizeof(double));
    if (new->x == NULL) {
        particle_destroy(new);
        return NULL;
    }

    new->velocity = (double *)malloc(nr_dimensions * sizeof(double));
    if (new->velocity == NULL) {
        particle_destroy(new);
        return NULL;
    }

    new->previous_best_x = (double *)malloc(nr_dimensions * sizeof(double));
    if (new->previous_best_x == NULL) {
        particle_destroy(new);
        return NULL;
    }

    new->neighbours = (particle_t **)malloc(nr_neighbours * sizeof(particle_t *));
    if (new->neighbours == NULL) {
        particle_destroy(new);
        return NULL;
    }

    new->neighbourhood_best_x = (double *)malloc(nr_dimensions * sizeof(double));
    if (new->neighbourhood_best_x == NULL) {
        particle_destroy(new);
        return NULL;
    }

    return new;
}

void particle_destroy(particle_t *particle)
{
    free(particle->neighbourhood_best_x);
    free(particle->neighbours);
    free(particle->previous_best_x);
    free(particle->velocity);
    free(particle->x);
    free(particle);
}

void particle_init(particle_t *particle, limits_t *dimensions_limits, double (*fitness_func)(double *x))
{
    int i;

    for (i = 0; i < dimensions_limits->nr_dimensions; ++i) {
        int min_d = dimensions_limits->dimensions[i].min;
        int max_d = dimensions_limits->dimensions[i].max;

        particle->x[i] = rd_uniform(min_d, max_d);
        particle->previous_best_x[i] = particle->x[i];
        particle->neighbourhood_best_x[i] = particle->x[i];

        particle->velocity[i] = rd_uniform(min_d - particle->x[i], max_d - particle->x[i]);
    }
    particle->fitness = fitness_func(particle->x);
    particle->previous_best_fitness = particle->fitness;
    particle->neighbourhood_best_fitness = particle->fitness;
}

bool set_neighbourhoods(particle_t **particles, pso_config_t *pso_config)
{
    int i, j, informed_idx;
    particle_t **best_neighbours;

    int nr_particles = pso_config->nr_particles;
    int k = pso_config->nr_neighbours;
    int nr_dimensions = pso_config->search_space_limits->nr_dimensions;

    best_neighbours = (particle_t**)malloc(nr_particles * sizeof(particle_t*));
    if (best_neighbours == NULL)
        return false;

    // at least informs itself
    for (i = 0; i < nr_particles; ++i)
        best_neighbours[i] = particles[i];

    // selects the best informants
    for (i = 0; i < nr_particles; ++i) {
        for (j = 0; j < k; ++j) {
            informed_idx = rd_iuniform(0, nr_particles);
            particles[i]->neighbours[j] = particles[informed_idx];
            if (particles[i]->previous_best_fitness < best_neighbours[informed_idx]->previous_best_fitness)
                best_neighbours[informed_idx] = particles[i];
        }
    }

    // copy the informants data
    // so the copy is done only once
    for (i = 0; i < nr_particles; ++i) {
        particles[i]->neighbourhood_best_fitness = best_neighbours[i]->previous_best_fitness;
        memcpy(particles[i]->neighbourhood_best_x, best_neighbours[i]->previous_best_x, nr_dimensions * sizeof(double));
    }

    free(best_neighbours);
    return true;
}

swarm_t *swarm_create(pso_config_t *pso_config, double (*fitness_func)(double *x))
{
    int i;
    swarm_t *new;
    particle_t *new_particle;

    new = (swarm_t *)calloc(1, sizeof(swarm_t));
    if (new == NULL)
        goto fail;

    new->particles = (particle_t **)calloc(pso_config->nr_particles, sizeof(particle_t *));
    if (new->particles == NULL)
        goto fail;

    // create and init each particle
    mt_seed();
    for (i = 0; i < pso_config->nr_particles; ++i) {
        new_particle = particle_create(pso_config);
        if (new_particle == NULL)
            goto fail;
        particle_init(new_particle, pso_config->search_space_limits, fitness_func);
        new->particles[i] = new_particle;
    }
    new->nr_particles = pso_config->nr_particles;

    new->search_space = pso_config->search_space_limits;

    if (set_neighbourhoods(new->particles, pso_config))
        return new;

fail:
    swarm_destroy(new);
    return NULL;
}

void swarm_destroy(swarm_t *swarm)
{
    int i;
    if (swarm->particles != NULL) {
        for (i = 0; i < swarm->nr_particles; ++i)
            particle_destroy(swarm->particles[i]);
    }
    free(swarm);
}

int iterarions(swarm_t *swarm, pso_config_t *pso_config, int nr_iterations)
{
    int i;
    double new_best_fitness;

    // update particles positions
    for (i = 0; i < nr_iterations; ++i) {
    }

    // update particles neighbourhoods
    if (new_best_fitness < swarm->best_fitness)
        set_neighbourhoods(swarm->particles, pso_config);
    else {
        // particles get informed by its neighbours
    }

    for (i = 0; i < nr_iterations; ++i) {
    }

    return i;
}
