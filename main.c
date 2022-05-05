#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "graph.h"
#include "util.h"

struct Graph read_in_file(char* filename);
void write_out_file(char* file_name, int* vector, int vector_size);

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
    int POPULATION_SIZE = 400; // 100, 200, 300, 400, ..., 1000
    double SELECTION_PRESSURE = 3.2; // x10 (3 ~ 4)
    double CROSSOVER_THRESHOLD = 0.90; // x10 (0 ~ 1)
    double EXECUTION_TIME = 177.0;

    // Init randomness
    srand(time(NULL));

    // Read graph data
    struct Graph graph_data = read_in_file(in_file_name);
    int** weight_table = graph_data.weight_table;
    int num_of_vertex = graph_data.num_of_vertex;

    // Init population data
    int** solutions = (int**)malloc(POPULATION_SIZE * sizeof(int*));
    double* fitnesses = (double*)malloc(POPULATION_SIZE * sizeof(double));
    int* values = (int*)malloc(POPULATION_SIZE * sizeof(int));
    for(int i = 0; i < POPULATION_SIZE; i++) {
        solutions[i] = (int*)malloc(num_of_vertex * sizeof(int));
        for (int j = 0; j < num_of_vertex; ++j)
            solutions[i][j] = 0;
        int iter = rand() % num_of_vertex;
        for (int j = 0; j < iter; ++j)
            solutions[i][rand() % num_of_vertex] = 1;
        // Generate the set of solutions with random values
//        for(int j = 0; j < num_of_vertex; j++)
//            solutions[i][j] = (int)(((double)rand() / ((double)(RAND_MAX) + 1.0)) * 2);
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
        int ci = values[i];
        int cw = worst_value;
        int cb = best_value;

        fitnesses[i] = (double)(ci - cw) + (cb - cw) / (SELECTION_PRESSURE - 1.0);
        sum_of_fitnesses += fitnesses[i];
    }

    clock_t start = clock();
    while(1) {
        clock_t now = clock();
        double time_spent = (double)(now - start) / CLOCKS_PER_SEC;

        if (time_spent > EXECUTION_TIME)
            break;

        // Selection
        int idx_of_mother = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);
        int idx_of_father = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);

        for (int i = 0; i < 10 && idx_of_mother == idx_of_father; ++i)
            idx_of_father = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);

        // Crossover
        int* child = (int*) malloc(sizeof (int) * num_of_vertex);
        for (int i = 0; i < num_of_vertex; ++i) {
            double r = rand() / ((double)RAND_MAX + 1.0);
            child[i] = r < CROSSOVER_THRESHOLD
                    ? solutions[idx_of_mother][i]
                    : solutions[idx_of_father][i];
        }

        // Mutation
        int mutated_index = (int) (rand() / ((double)RAND_MAX + 1.0) * num_of_vertex) % num_of_vertex;
        child[mutated_index] = !child[mutated_index];

        // Replace with worst case
        for (int i = 0; i < num_of_vertex; ++i)
            solutions[worst_solution_index][i] = child[i];

        // Before update value, reduce sum of fitness
        sum_of_fitnesses -= fitnesses[worst_solution_index];

        // Update value of replaced solution
        values[worst_solution_index] = evaluate(graph_data, child);
        fitnesses[worst_solution_index] = (double)(values[worst_solution_index] - worst_value) + (best_value - worst_value) / (SELECTION_PRESSURE - 1.0);

        // Update sum of fitness
        sum_of_fitnesses += fitnesses[worst_solution_index];

        SAFE_FREE(child);

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

    for(int i = 0; i < POPULATION_SIZE; i++)
        free(solutions[i]);
    SAFE_FREE(solutions);

    for(int i = 0; i < num_of_vertex; i++)
        free(weight_table[i]);
    SAFE_FREE(weight_table);

    return 0;
}

struct Graph read_in_file(char* filename) {
    int num_of_vertex, num_of_edges;
    int **graph = NULL;
    FILE* in_file = fopen(filename, "r");
    if (in_file == NULL) {
        printf("File not found %s\n", filename);
        exit(-1);
    }

    fscanf(in_file, "%d %d", &num_of_vertex, &num_of_edges);

    // Initialize graph
    graph = (int**) malloc(sizeof(int*) * num_of_vertex);
    for (int i = 0; i < num_of_vertex; ++i) {
        graph[i] = (int *) malloc(sizeof(int) * num_of_vertex);
        for (int j = 0; j < num_of_vertex; ++j) {
            graph[i][j] = 0;
        }
    }

    // Mark weight of graph
    for(int i = 0; i < num_of_edges; i++) {
        int from, to, value;
        fscanf(in_file, "%d %d %d", &from, &to, &value);
        graph[from-1][to-1] = value;
        graph[to-1][from-1] = value;
    }

    fclose(in_file);

    return init_graph(num_of_vertex, num_of_edges, graph);
}

void write_out_file(char* file_name, int* vector, int vector_size) {
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
