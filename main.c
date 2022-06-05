#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "graph.h"
#include "util.h"

struct Graph read_in_file(char* filename);
void write_out_file(char* file_name, const int* vector, int vector_size);

/**
 * Main func
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {

    char* in_file_name = argv[1];
    char* out_file_name = argv[2];

    // Hyper parameters
    int POPULATION_SIZE = 100; // 100, 200, 300, 400, ..., 1000
    double SELECTION_PRESSURE = 3.0; // x10 (3 ~ 4)

    // Init randomness
    srand(time(NULL));

    // Read graph data
    struct Graph graph_data = read_in_file(in_file_name);
    int** edges = graph_data.edges;
    int num_of_vertex = graph_data.num_of_vertex;

    double EXECUTION_TIME = 175.0;
    EXECUTION_TIME = num_of_vertex / 6 - 5;

    // Init population data
    int** solutions = (int**)malloc(POPULATION_SIZE * sizeof(int*));
    double* fitnesses = (double*)malloc(POPULATION_SIZE * sizeof(double));
    int* values = (int*)malloc(POPULATION_SIZE * sizeof(int));
    for(int i = 0; i < POPULATION_SIZE; ++i) {
        solutions[i] = (int*)malloc(num_of_vertex * sizeof(int));
        // Generate the set of solutions with random values
        for(int j = 0; j < num_of_vertex; ++j)
            solutions[i][j] = rand() % 2;
        fitnesses[i] = 0;
        values[i] = evaluate(graph_data, solutions[i]);
    }

    // update worst and best values
    struct MinAvgMax min_avg_max = get_min_avg_max_from_vector(values, POPULATION_SIZE);

    int worst_solution_index = min_avg_max.min_idx;
    int best_solution_index = min_avg_max.max_idx;

    int worst_value = values[worst_solution_index];
    int best_value = values[best_solution_index];

    // Update fitness of each solution
    double sum_of_fitnesses = 0;
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        fitnesses[i] = (double)(values[i] - worst_value) + (best_value - worst_value) / (SELECTION_PRESSURE - 1.0);
        sum_of_fitnesses += fitnesses[i];
    }

    clock_t start = clock();
    while(1) {
        clock_t now = clock();
        double time_spent = (double)(now - start) / CLOCKS_PER_SEC;

        if (time_spent > EXECUTION_TIME)
            break;

        // Do local optimization
        int temp=rand()%POPULATION_SIZE;
        local_opt(graph_data, solutions[temp]);

        // Update value
        values[temp] = evaluate(graph_data, solutions[temp]);

        // Update best, worst case solution
        min_avg_max = get_min_avg_max_from_vector(values, POPULATION_SIZE);

        worst_solution_index = min_avg_max.min_idx;
        best_solution_index = min_avg_max.max_idx;
        best_value = values[best_solution_index];

       double avg_of_values = min_avg_max.avg_value;
       printf("%.2f, %d, %.2f\n", time_spent, best_value, avg_of_values);
    }

    // Write output file
    write_out_file(out_file_name, solutions[best_solution_index], num_of_vertex);

    // Free allocated memory
    for(int i = 0; i < POPULATION_SIZE; i++)
        MACRO_FREE(solutions[i]);
    MACRO_FREE(solutions);

    MACRO_FREE(values);
    MACRO_FREE(fitnesses);

    for(int i = 0; i < num_of_vertex; i++)
        MACRO_FREE(edges[i]);
    MACRO_FREE(edges);

    return 0;
}

struct Graph read_in_file(char* filename) {
    FILE* in_file = fopen(filename, "r");
    if (in_file == NULL) {
        printf("File not found %s\n", filename);
        exit(-1);
    }

    int num_of_vertex, num_of_edges;
    fscanf(in_file, "%d %d", &num_of_vertex, &num_of_edges);

    // Initialize graph
    int **graph = (int**) malloc(sizeof(int*) * num_of_vertex);
    for (int i = 0; i < num_of_vertex; ++i) {
        graph[i] = (int *) malloc(sizeof(int) * num_of_vertex);
        for (int j = 0; j < num_of_vertex; ++j)
            graph[i][j] = 0;
    }

    // Mark weight of graph
    int from, to, value;
    for(int i = 0; i < num_of_edges; ++i) {
        fscanf(in_file, "%d %d %d", &from, &to, &value);
        graph[from-1][to-1] = graph[to-1][from-1] =value;
    }

    fclose(in_file);

    return init_graph(num_of_vertex, num_of_edges, graph);
}

void write_out_file(char* file_name, const int* vector, int vector_size) {
    // Write output file
    FILE* out_file = fopen(file_name, "w");
    if (out_file == NULL) {
        printf("Something wrong while opening output file.\n");
        exit(-1);
    }

    int flag = vector[0];
    for (int i = 0; i < vector_size; ++i)
        // Warning: start index from '1'
        if (vector[i] == flag)
            fprintf(out_file, "%d ", i + 1);
    fclose(out_file);
}
