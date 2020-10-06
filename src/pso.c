#include "pso.h"

#include <stdlib.h>

#include <randistrs.h>

particle_t *particle_create(limits_t *dimensions_limits)
{
    particle_t *new;
    int nr_dimensions = dimensions_limits->nr_dimensions;

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

    return new;
}

void particle_destroy(particle_t *particle)
{
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
        particle->previous_best[i] = particle->x[i];
        particle->velocity[i] = rd_uniform(min_d - particle->x[i], max_d - particle->x[i]);
    }
    particle->fitness = fitness_func(particle->x);
    particle->previous_best_fitness = particle->fitness;
    particle->neighbourhood_best_fitness = particle->fitness;
}

void set_neighbourhoods(particle_t **particles, int nr_particles, int k)
{
    int i, j;
    particle_t *particle_informed;

    for (i = 0; i < nr_particles; ++i) {
        // at least informs itself
        particles[i]->neighbourhood_best_fitness = particles[i]->fitness;
        for (j = 0; j < k; ++j) {
            particle_informed = particles[rd_iuniform(0, nr_particles)];
            particles[i]->neighbours[j] = particle_informed;
            if (particle_informed->fitness < particles[i]->neighbourhood_best_fitness)
                particles[i]->neighbourhood_best_fitness = particle_informed->fitness;
        }
    }
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
        new_particle = particle_create(pso_config->search_space_limits);
        if (new_particle == NULL)
            goto fail;
        particle_init(new_particle, pso_config->search_space_limits, fitness_func);
        new->particles[i] = new_particle;
    }
    new->nr_particles = pso_config->nr_particles;

    new->search_space = pso_config->search_space_limits;

    set_neighbourhoods(new->particles, new->nr_particles, pso_config->nr_neighbours);

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

int iterarions(swarm_t *swarm, int nr_iterations)
{
    int i;

    // update particles positions
    for (i = 0; i < nr_iterations; ++i) {
    }

    // update particles neighbourhoods
    for (i = 0; i < nr_iterations; ++i) {
    }

    return i;
}
