#include "pso.h"

#include <stdlib.h>

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
