#include <parallel_evolution.h>
#include <math.h>
#include <stdlib.h>
#include "pso.h"

#define ERR_MALLOC 1

#define MODULE_APP "pe-pso"

#define A 10    /* arbitrary constant, default for rastrigin */
/*#define N 50    /* number of dimensions */

/* define here the fitness function that will be used by all algorithms of parallel_evolution */
/* fitness function for generalized rastrigin function */
double objective_function(double *x)
{
    double sum = 0;
    int i, n;

    n = parallel_evolution_get_number_of_dimensions();

    for (i = 0; i < n; ++i)
        sum += x[i] * x[i] - A * cos(2 * M_PI * x[i]);

    return A * n + sum;
}

limits_t limits;
pso_config_t pso_config;
swarm_t *swarm;
int pso_iterations = 0;
int pso_max_iterations = 5000;  // default value

void pso_init(config_t *config)
{
    int i, n;
    dimension_t *dimensions = NULL;

    parallel_evolution_config_lookup_int(config, "pso.max_iterations", &pso_max_iterations);

    n = parallel_evolution_get_number_of_dimensions();
    dimensions = (dimension_t *)malloc(n * sizeof(dimension_t));    // TODO: dealloc this somewhere
    if (dimensions == NULL) {
        parallel_evolution_log(LOG_PRIORITY_ERR, MODULE_APP, "Couldn't allocate memory for dimensions limits. This is the end...");
        exit(ERR_MALLOC);
    }
    for (i = 0; i < n; ++i) {
        dimensions[i].min = parallel_evolution_get_limit_min(i);
        dimensions[i].max = parallel_evolution_get_limit_max(i);
    }
    limits.dimensions = dimensions;
    limits.nr_dimensions = n;

    parallel_evolution_config_lookup_int(config, "pso.nr_particles", &pso_config.nr_particles);
    parallel_evolution_config_lookup_int(config, "pso.nr_neighbours", &pso_config.nr_neighbours);
    pso_config.search_space_limits = &limits;

    swarm = swarm_create(&pso_config, objective_function);
}

void pso_run_iterations(int nr_iterations)
{
    pso_iterations += iterations(swarm, &pso_config, objective_function, nr_iterations);
}

void pso_insert_migrant(migrant_t *migrant)
{
    int i, worse_i, d, k;
    particle_t *worse_particle;

    // find the particle with the worse fitness
    for (worse_i = 0, i = 1; i < pso_config.nr_particles; ++i)
        if (swarm->particles[i]->fitness > swarm->particles[worse_i]->fitness)
            worse_i = i;
    worse_particle = swarm->particles[worse_i];

    // replace the position of that particle
    for (d = 0; d < limits.nr_dimensions; ++d)
        worse_particle->x[d] = migrant->var[d];

    // compute fitness
    worse_particle->fitness = objective_function(worse_particle->x);

    // update the previous best and its neighbourhood
    if (worse_particle->fitness < worse_particle->previous_best_fitness) {
        worse_particle->previous_best_fitness = worse_particle->fitness;
        for (k = 0; k < pso_config.nr_neighbours; ++k) {
            if (worse_particle->previous_best_fitness < worse_particle->neighbours[k]->neighbourhood_best_fitness) {
                worse_particle->neighbours[k]->neighbourhood_best_fitness = worse_particle->previous_best_fitness;
                worse_particle->neighbours[k]->neighbourhood_best_x = worse_particle->previous_best_x;
            }
        }
    }
}

void pso_pick_migrant(migrant_t *migrant)
{
    int i, best_i;
    particle_t *best_particle;

    for (best_i = 0, i = 1; i < pso_config.nr_particles; ++i)
        if (swarm->particles[i]->fitness < swarm->particles[best_i]->fitness)
            best_i = i;
    best_particle = swarm->particles[best_i];

    for (i = 0; i < migrant->var_size; ++i)
        migrant->var[i] = best_particle->x[i];
}

int pso_ended()
{
    return pso_iterations >= pso_max_iterations;
}

algorithm_stats_t *pso_get_stats()
{
    int i;
    algorithm_stats_t *stats;
    double avg_fitness = 0;

    stats = (algorithm_stats_t *)malloc(sizeof(algorithm_stats_t));
    if (stats == NULL)
        return NULL;

    for (i = 0; i < swarm->nr_particles; ++i) {
        avg_fitness += swarm->particles[i]->fitness;
    }
    avg_fitness /= swarm->nr_particles;

    stats->iterations = swarm->iteration;
    stats->fitness_evals = swarm->iteration * swarm->nr_particles;
    stats->best_fitness = swarm->best_fitness;
    stats->avg_fitness = avg_fitness;

    return stats;
}

status_t pso_get_population(population_t **pop2send)
{
    int i;

    if (population_create(pop2send, swarm->nr_particles) != SUCCESS) {
        return FAIL;
    }

    for (i = 0; i < swarm->nr_particles; ++i) {
        (*pop2send)->individuals[i]->var = swarm->particles[i]->x;
        (*pop2send)->individuals[i]->var_size = swarm->search_space->nr_dimensions;
    }

    return SUCCESS;
}

int main(int argc, char **argv)
{
    algorithm_t *pso;
    int ret;
    algorithm_create(&pso,
                        pso_init,
                        pso_run_iterations,
                        pso_insert_migrant,
                        pso_pick_migrant,
                        pso_ended,
                        pso_get_population,   // the population that will be sent to the main node
                        pso_get_stats);
    parallel_evolution_set_algorithm(pso);
    ret = parallel_evolution_run(&argc, &argv);

    algorithm_destroy(&pso);

    return ret;
}
