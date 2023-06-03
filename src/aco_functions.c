#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int pick_next_city(int *visited_cities, float **dist, float **pheromone, int pheromonePower, int distPower, int n_cities, int id){
    // This function picks the next city to visit based on the pheromone trails and the distance to the cities
    float prob[n_cities];
    float temp;
    float sum = 0;
    for(int i=0; i<n_cities; i++){
        if (visited_cities[i] == 1) // Only consider unvisited cities
        {
            prob[i] = 0;
        }
        else{
            temp = pow(pheromone[id][i],pheromonePower) * pow(1/dist[id][i],distPower);
            sum += temp;
            prob[i] = sum;
        }
    }
    // Normalize probabilities
    for(int i=0; i<n_cities; i++){
        if (visited_cities[i] == 0)
        {
            prob[i] = prob[i]/sum;
        }
    }
    // Pick a city based on probabilities
    float u = (float)rand()/(float)(RAND_MAX); // Random number between 0. and 1.
    int j = 0;
    while (u > prob[j])
    {
        j++;
    }
    return j;
}

struct Tour{
    float length;
    int* path;
    
};

struct Tour let_ants_do_a_tour(int L, int n_cities, float **dist, float **pheromone, int pheromonePower, int distPower, struct Tour best_tour){  
        // This function contructs tours for L ants
        int j, i;
        for(int a = 0; a<L; a++){
            // Reset tour and unvisited cities
            struct Tour tour = {0., (int*) calloc(n_cities+1,sizeof(int))};
            int *visited_cities = (int*) calloc(n_cities,sizeof(int));
            int id = 0;
            // Start in City 0
            visited_cities[0] = 1;
            tour.path[0] = 0;
            // Build the path
            for(j=1; j<n_cities; j++){
                // Select next city by determining probabilities and doing a MC selection
                int id_new = pick_next_city(visited_cities, dist, pheromone, pheromonePower, distPower, n_cities, id);
                // Remove the selected city from list of unvisited cities
                visited_cities[id_new] = 1;
                // Update tour
                tour.path[j] = id_new;
                tour.length += dist[id][id_new];
                id = id_new;
            }
            // Return to the origin
            tour.path[n_cities] = 0;
            tour.length += dist[id][0];
            // Update best solution found in this iteration
            if(tour.length < best_tour.length){
                best_tour.length = tour.length;
                for(i=0; i<n_cities; i++){
                    best_tour.path[i] = tour.path[i];
                }
            }
            free(tour.path);
            free(visited_cities);
        }

    return best_tour;
}

float** update_pheromone(int n_cities, struct Tour best_tour, float **pheromone, float evaporationRate, float pheromone_min, float pheromone_max){
    // This function updates the pheromone matrix
    int i, j;

    // Initialize delta_pheromone
    float **delta_pheromone = (float **) malloc(n_cities * sizeof(float *));
    delta_pheromone[0] = (float *) malloc(n_cities * n_cities * sizeof(float));
    for ( i = 1; i < n_cities; i++) {
    delta_pheromone[i] = delta_pheromone[0] + i * n_cities;
    }

    // Compute delta_pheromone based on the best tour
    for(i=0; i<n_cities-1; i++){
        delta_pheromone[best_tour.path[i]][best_tour.path[i+1]] = 1/best_tour.length;
    }

    // Update pheromone matrix
    for(i=0; i<n_cities; i++) {
        for(j=0; j<n_cities; j++) {
            float temp = (1-evaporationRate) * pheromone[i][j] + delta_pheromone[i][j];
            // Keep within bounds
            if(temp > pheromone_max){
                pheromone[i][j] = pheromone_max;
            }
            else if(temp < pheromone_min){
                pheromone[i][j] = pheromone_min;
            }
            else{
                pheromone[i][j] = temp;
            }
        }
    }
    return pheromone;
}

void write_tour_to_file(const char * file_path, char w, int n_cities, struct Tour tour_to_export){
    int i;
    FILE *fp;

    char delim[] = "\0";
    const char *file_path_tmp = strtok((char *)file_path, delim);
    
    fp = fopen(file_path_tmp, &w);
    for (i = 0; i < n_cities; i++) {
        fprintf(fp, "%d, ", tour_to_export.path[i]);
    }
    fprintf(fp, "%f\n", tour_to_export.length);

    fclose(fp);
}

void write_pheromone_to_file(const char * file_path, char w, int n_cities, float **pheromone){
    int i, j;
    FILE *fp;

    char delim[] = "\0";
    const char *file_path_tmp = strtok((char *)file_path, delim);
    
    fp = fopen(file_path_tmp, &w);
    for (i = 0; i < n_cities; i++) {
        for (j = 0; j < n_cities; j++){
            if(i == n_cities-1 && j == n_cities-1){
                fprintf(fp, "%f", pheromone[i][j]);
            }else{
            fprintf(fp, "%f, ", pheromone[i][j]);
            }            
        }
    }
    fprintf(fp, "\n");

    fclose(fp);
}

struct Params{
    int n_cities;
    int n_iterations;
    int n_ants;
    float evaporationRate;
    int pheromonePower;
    int distPower;
    float length_guess;
};

struct Params read_parameters(const char * param_file){
    FILE *fptr;
    char input_line[100];
    // Open the parameter file in read mode
    fptr = fopen(param_file, "r");
    struct Params params;
    if(fptr != NULL) {
        // Read the content of the file
        fgets(input_line, 100, fptr);
        char delim[] = "\0";
        char *ptr = strtok(input_line, delim);
        params.n_cities = atoi(ptr);
        fgets(input_line, 100, fptr);
        ptr = strtok(input_line, delim);
        params.n_iterations = atoi(ptr);
        fgets(input_line, 100, fptr);
        ptr = strtok(input_line, delim);
        params.n_ants = atoi(ptr);
        fgets(input_line, 100, fptr);
        ptr = strtok(input_line, delim);
        params.evaporationRate = atof(ptr);
        fgets(input_line, 100, fptr);
        ptr = strtok(input_line, delim);
        params.pheromonePower = atoi(ptr);
        fgets(input_line, 100, fptr);
        ptr = strtok(input_line, delim);
        params.distPower = atoi(ptr);
        fgets(input_line, 100, fptr);
        ptr = strtok(input_line, delim);
        params.length_guess = atof(ptr);
        fclose(fptr);
    } else {
        printf("Not able to open parameter file.");
    }
    return params;
}

float** read_locations(const char * location_file, int n_cities){
    int i;
    float **locations = (float **) malloc(n_cities * sizeof(float *));
    locations[0] = (float *) malloc(n_cities * 2 * sizeof(float));
    for (i = 1; i < n_cities; i++) {
        locations[i] = locations[0] + i * 2;
    }
    FILE *fptr;
    char input_line[100];
    // Open the location file in read mode
    fptr = fopen(location_file, "r");
    // If the file exist
    if(fptr != NULL) {
        // Read the content of the file
        for(i=0; i<n_cities; i++){
            fgets(input_line, 100, fptr);
            char delim[] = ",";
            char *ptr = strtok(input_line, delim);
            locations[i][0] = atof(ptr);
            ptr = strtok(NULL, delim);
            locations[i][1] = atof(ptr);
        }
        fclose(fptr);
    } else {
        printf("Not able to open parameter file.");
    }
    return locations;
}