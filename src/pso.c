#include "pso.h"

#include <stdlib.h>

#define MT_NO_INLINE
#include <mtwist.h>

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

double rnd_in_range(double min, double max)
{
    return min + mt_ldrand() * (max - min);
}

void particle_init(particle_t *particle, limits_t *dimensions_limits, double (*fitness_func)(double *x))
{
    for (int i = 0; i < dimensions_limits->nr_dimensions; ++i) {
        int min_d = dimensions_limits->dimensions[i].min;
        int max_d = dimensions_limits->dimensions[i].max;

        particle->x[i] = rnd_in_range(min_d, max_d);
        particle->previous_best[i] = particle->x[i];
        particle->velocity[i] = rnd_in_range(min_d - particle->x[i], max_d - particle->x[i]);
    }
    particle->fitness = fitness_func(particle->x);
    particle->previous_best_fitness = particle->fitness;
}

swarm_t *swarm_create(limits_t *dimensions_limits, int nr_particles, double (*fitness_func)(double *x))
{
    swarm_t *new;
    particle_t *new_particle;

    new = (swarm_t *)calloc(1, sizeof(swarm_t));
    if (new == NULL)
        goto fail;

    new->particles = (particle_t **)calloc(nr_particles, sizeof(particle_t *));
    if (new->particles == NULL)
        goto fail;
    new->nr_particles = nr_particles;

    // create and init each particle
    mt_seed();
    for (int i = 0; i < nr_particles; ++i) {
        new_particle = particle_create(dimensions_limits);
        if (new_particle == NULL)
            goto fail;
        particle_init(new_particle, dimensions_limits, fitness_func);
    }

    new->search_space = dimensions_limits;

    return new;

fail:
    swarm_destroy(new);
    return NULL;
}

void swarm_destroy(swarm_t *swarm)
{
    if (swarm->particles != NULL) {
        for (int i = 0; i < swarm->nr_particles; ++i)
            particle_destroy(swarm->particles[i]);
    }
    free(swarm);
}
