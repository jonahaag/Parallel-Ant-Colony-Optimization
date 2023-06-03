// Ant colony optimization
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "aco_functions.c"

int main(int argc, char const *argv[])
{

    int tag, P, p, rc, step;
    tag = 42;

    MPI_Status status;
    rc = MPI_Init(NULL, NULL);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &P);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &p);

    // Starting time
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    // Initialize TSP
    srand((unsigned) p);
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
    
    // Initialize locations
    float **locations = read_locations(location_file, n_cities);

    // Initialize best tour, only used on process with rank 0
    struct Tour best_tour_global = {1000000000., (int*) calloc(n_cities+1,sizeof(int))};

    // Initialize distance matrix and pheromone trails
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

    // Compute distance matrix and set initial pheromone trails
    for(i=0; i<n_cities; i++) {
        for(j=0; j<n_cities; j++) {
            dist[i][j] = sqrt(pow(locations[i][0]-locations[j][0],2) + pow(locations[i][1]-locations[j][1],2));
            pheromone[i][j] = pher_init;
        }
    }
    free(locations);

    //Distribute ants on processes 
    int L = n_ants/P;
    int R = n_ants%P;
    if(p < R){
        L +=1;
    }

    // For number of iterations or some stopping criterion
    for(it=0; it<n_iterations; it++){
            
        // Construct ant paths
        struct Tour best_tour_local = {1000000000., (int*) calloc(n_cities+1,sizeof(int))};

        best_tour_local = let_ants_do_a_tour(L, n_cities, dist, pheromone, pheromonePower, distPower, best_tour_local);

        // Send the best past length of each process to the master process
        float* rbuf = (float *)malloc(P*sizeof(float));
        MPI_Gather(&best_tour_local.length, 1, MPI_FLOAT, rbuf, 1, MPI_FLOAT, 0, MPI_COMM_WORLD); 

        // Check if local path is better than global path
        int min_id = -1;
        if (p == 0){
            for (i =0; i < P; i++){
                if (rbuf[i] < best_tour_global.length){
                    best_tour_global.length = rbuf[i];
                    min_id = i;
                }
            }
        }
        free(rbuf);

        // Get the path corresponding to the best (shortest) length
        MPI_Barrier(MPI_COMM_WORLD);
        rc = MPI_Bcast(&min_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (min_id == 0){ // The best path is on the process with rank 0
            for (i =0; i < n_cities+1; i++){
                best_tour_global.path[i] = best_tour_local.path[i];
            }
        } 
        else if (min_id > 0){ // The best path is on another process
            if (p == min_id) {
                // Send the best path to the master process 
                rc = MPI_Send(best_tour_local.path, n_cities+1, MPI_INT, 0, tag, MPI_COMM_WORLD); 
            }

            if (p == 0){  
                // Receive the best path
                rc = MPI_Recv(best_tour_global.path, n_cities+1, MPI_INT, min_id, tag, MPI_COMM_WORLD, &status);
            }
        }

        if (p==0){
            // Update pheromone trails based on solution found
            pheromone = update_pheromone(n_cities, best_tour_global, pheromone, evaporationRate, pheromone_min, pheromone_max);
        }

        // Send the updated pheromone matrix to all processes
        for(i=0; i<n_cities; i++) {
            rc = MPI_Bcast(pheromone[i], n_cities, MPI_FLOAT, 0, MPI_COMM_WORLD);
        }
    }
    if (p==0){
        printf("Best path found: %f.\n", best_tour_global.length);
    }

    // ending time
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    if (p == 0){
        printf("Time elapsed: %f seconds", end-start);
    }
        
    MPI_Finalize();
    return 0;
}
