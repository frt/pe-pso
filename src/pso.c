#include "pso.h"

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

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

void inform_neighbours(particle_t **particles, pso_config_t *pso_config)
{
    int i, j;

    int nr_particles = pso_config->nr_particles;
    int k = pso_config->nr_neighbours;
    int nr_dimensions = pso_config->search_space_limits->nr_dimensions;

    // at least informs itself
    for (i = 0; i < nr_particles; ++i) {
        if (particles[i]->previous_best_fitness < particles[i]->neighbourhood_best_fitness) {
            particles[i]->neighbourhood_best_fitness = particles[i]->previous_best_fitness;
            memcpy(particles[i]->neighbourhood_best_x, particles[i]->previous_best_x, nr_dimensions * sizeof(double));
        }
    }

    // inform its neighbours
    for (i = 0; i < nr_particles; ++i) {
        for (j = 0; j < k; ++j) {
            if (particles[i]->previous_best_fitness < particles[i]->neighbours[j]->neighbourhood_best_fitness) {
                particles[i]->neighbours[j]->neighbourhood_best_fitness = particles[i]->previous_best_fitness;
                memcpy(particles[i]->neighbours[j]->neighbourhood_best_x, particles[i]->previous_best_x, nr_dimensions * sizeof(double));
            }
        }
    }
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

        if (i == 0 || new_particle->fitness < new->best_fitness)
            new->best_fitness = new_particle->fitness;
    }
    new->nr_particles = pso_config->nr_particles;
    new->best_fitness_iteration = 0;

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

// updates velocity and move particle
bool move_particle(particle_t *particle, int nr_dimensions, double c, double w)
{
    int i;
    double *gravity_center, *x1;
    double *x = particle->x;
    double *p = particle->previous_best_x;
    double *l = particle->neighbourhood_best_x;
    double radius, r, length;

    gravity_center = (double *)malloc(nr_dimensions * sizeof(double));
    x1             = (double *)malloc(nr_dimensions * sizeof(double));
    if (gravity_center == NULL || x1 == NULL) {
        free(gravity_center);
        free(x1);
        return false;
    }

    // calculate the gravity center
    for (i = 0; i < nr_dimensions; ++i) {
        gravity_center[i] = x[i] + c * ((p[i] + l[i] - 2 * x[i]) / 3);
    }

    // radius ||Gi − xi|| around gravity_center makes a hypersphere
    for (i = 0, radius = 0; i < nr_dimensions; ++i)
        radius += pow((gravity_center[i] - x[i]), 2);
    radius = sqrt(radius);

    /* get a random point x1 inside the hypersphere
       [https://baezortega.github.io/2018/10/14/hypersphere-sampling/]
    */
    // 1. define the direction
    for (i = 0, length = 0; i < nr_dimensions; ++i) {
        x1[i] = rd_normal(0, 1);    // μ = 0 and σ = 1
        length += pow(x1[i], 2);
    }
    length = sqrt(length);

    // 2. define radius normalized in [0, 1)
    r = pow(rd_uniform(0, 1), 1 / nr_dimensions);

    // 3. apply the random radius and move the point toward the gravity center
    for (i = 0; i < nr_dimensions; ++i)
        x1[i] = gravity_center[i] + radius * r * (x1[i] / length);

    // vi(t+1) = w * vi(t) + x'i(t) − xi(t)
    // xi(t+1) = xi(t) + vi(t+1)
    for (i = 0; i < nr_dimensions; ++i) {
        particle->velocity[i] = w * particle->velocity[i] + x1[i] - x[i];
        particle->x[i] += particle->velocity[i];
    }

    free(gravity_center);
    free(x1);
    return true;
}

int iterarions(swarm_t *swarm, pso_config_t *pso_config, double (*fitness_func)(double *x), int nr_iterations)
{
    int i, j;
    double new_best_fitness = swarm->best_fitness;
    int nr_dimensions = swarm->search_space->nr_dimensions;
    int nr_particles = swarm->nr_particles;

    /* Maurice Clerc.
       Stagnation analysis in particle swarm optimization or
       what happens when nothing happens, http://hal.archives-ouvertes.fr/hal-00122031.
       Technical report, 2006. */
    static const double w = 1 / (2 * log(2));
    static const double c = 1 / 2 + log(2);

    // update particles positions
    for (j = 0; j < nr_iterations; ++j) {
        for (i = 0; i < nr_particles; ++i) {
            move_particle(swarm->particles[i], nr_dimensions, c, w);

            // calculate new fitness, etc.
            swarm->particles[i]->fitness = fitness_func(swarm->particles[i]->x);
            if (swarm->particles[i]->fitness < swarm->particles[i]->previous_best_fitness) {
                swarm->particles[i]->previous_best_fitness = swarm->particles[i]->fitness;
                memcpy(swarm->particles[i]->previous_best_x, swarm->particles[i]->x, nr_dimensions * sizeof(double));

                if (swarm->particles[i]->fitness < new_best_fitness)
                    new_best_fitness = swarm->particles[i]->fitness;
            }
        }

        // update particles neighbourhoods
        if (new_best_fitness < swarm->best_fitness) {
            swarm->best_fitness = new_best_fitness;
            inform_neighbours(swarm->particles, pso_config);    // particles get informed by its neighbours
        } else {
            if (!set_neighbourhoods(swarm->particles, pso_config))  // some error occurred, stop and return the number of iterations run until now
                return j;
        }
    }

    return j;
}
