// Ant colony optimization
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "aco_functions.c"

int main(int argc, char const *argv[])
{

    clock_t start, end;
    double time_taken;
    start = clock();

    // Initialize TSP
    int read_cities = 0;
    int i, j, it;

    const char * location_file = argv[1];
    const char * param_file = argv[2];

    // Read parameters from file
    struct Params params = read_parameters(param_file);
    int n_cities = params.n_cities;
    int n_iterations = params.n_iterations;
    int n_ants = params.n_ants;
    float evaporationRate = params.evaporationRate;
    int pheromonePower = params.pheromonePower;
    int distPower = params.distPower;
    float pheromone_max = 1./ (evaporationRate * params.length_guess);
    float pheromone_min = pheromone_max / (2. * (float) n_cities);
    float pher_init = pheromone_max;
    
    // Initialize locations, distance, and pheromone matrix
    float **locations = read_locations(location_file, n_cities);

    struct Tour best_tour = {1000000000, (int*) calloc(n_cities+1,sizeof(int))};

    float **dist = (float **) malloc(n_cities * sizeof(float *));
    dist[0] = (float *) malloc(n_cities * n_cities * sizeof(float));
    for ( i = 1; i < n_cities; i++) {
        dist[i] = dist[0] + i * n_cities;
    }

    float **pheromone = (float **) malloc(n_cities * sizeof(float *));
    pheromone[0] = (float *) malloc(n_cities * n_cities * sizeof(float));
    for ( i = 1; i < n_cities; i++) {
        pheromone[i] = pheromone[0] + i * n_cities;
    }

    // Compute distance matrix and initialize pheromone trails
    for(i=0; i<n_cities; i++) {
        for(j=0; j<n_cities; j++) {
            dist[i][j] = sqrt(pow(locations[i][0]-locations[j][0],2) + pow(locations[i][1]-locations[j][1],2));
            pheromone[i][j] = pher_init;
        }
    }
    free(locations);

    // For number of iterations or some stopping criterion
    for(it=0; it<n_iterations; it++){

        // Contruct ant paths and directly update the best tour
        best_tour = let_ants_do_a_tour(n_ants, n_cities, dist, pheromone, pheromonePower, distPower, best_tour);

        // Update pheromone trails based on solution found
        pheromone = update_pheromone(n_cities, best_tour, pheromone, evaporationRate, pheromone_min, pheromone_max);

    }
    printf("Best path found: %f.\n", best_tour.length);

    end = clock();
    time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Computation time = %lf seconds", time_taken);
    return 0;
}
