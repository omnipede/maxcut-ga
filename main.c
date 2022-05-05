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
    int POPULATION_SIZE = 400;
    double CROSSOVER_THRESHOLD = 0.9;
    int OFFSPRING_SIZE = 50;

    double EXECUTION_TIME = 177.0;

    // Init randomness
    srand(time(NULL));

    // Read graph data
    struct Graph graph_data = read_in_file(in_file_name);
    int** weight_table = graph_data.weight_table;
    int num_of_vertex = graph_data.num_of_vertex;

    // Init population data
    int** solutions = (int**)malloc(POPULATION_SIZE * sizeof(int*));
    int* values = (int*)malloc(POPULATION_SIZE * sizeof(int));
    for(int i = 0; i < POPULATION_SIZE; i++) {
        solutions[i] = (int*)malloc(num_of_vertex * sizeof(int));
        // Generate the set of solutions with random values
        for(int j = 0; j < num_of_vertex; j++)
            solutions[i][j] = (int)(((double)rand() / ((double)(RAND_MAX) + 1.0)) * 2);
        values[i] = evaluate(graph_data, solutions[i]);
    }

    // Init offspring buffer
    int** offsprings = (int**) malloc( POPULATION_SIZE * sizeof(int*));
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        offsprings[i] = (int*) malloc(num_of_vertex * sizeof(int));
        for (int j = 0; j < num_of_vertex; ++j)
            offsprings[i][j] = 0;
    }

    // update worst and best values
    struct MinAvgMax min_avg_max = get_min_avg_max_from_vector(values, POPULATION_SIZE);

    int best_solution_index = min_avg_max.max_idx;
    int best_value = values[best_solution_index];

    // Update fitness of each solution
    int generation = 0;
    clock_t start = clock();
    while(1) {
        generation += 1;
        if (generation > 2000) {
            OFFSPRING_SIZE = POPULATION_SIZE;
            // TODO selection count ?
        }

        clock_t now = clock();
        double time_spent = (double)(now - start) / CLOCKS_PER_SEC;

        if (time_spent > EXECUTION_TIME)
            break;

        // Selection
        int* parents = (int*) malloc(sizeof (int) * 2);
        select_two_from_vector(values, POPULATION_SIZE, parents);

        int idx_of_mother = parents[0];
        int idx_of_father = parents[1];
        SAFE_FREE(parents);

        // Crossover
        for (int i = 0; i < OFFSPRING_SIZE; ++i) {
            for (int j = 0; j < num_of_vertex; ++j) {
                double r = rand() / ((double)RAND_MAX + 1.0);
                offsprings[i][j] = r < CROSSOVER_THRESHOLD
                                   ? solutions[idx_of_mother][j]
                                   : solutions[idx_of_father][j];
            }
        }

        // Mutation
        if (generation > 200) {
            for (int i = 0; i < OFFSPRING_SIZE; ++i) {
                int mut_idx = rand() % num_of_vertex;
                offsprings[i][mut_idx] = !offsprings[i][mut_idx];
            }
        }

        // Replace with random
        for (int i = 0; i < OFFSPRING_SIZE; ++i) {
            // Take one from population to replace
            int replaced_idx = rand() % POPULATION_SIZE;
            while(replaced_idx == best_solution_index)
                replaced_idx = rand() % POPULATION_SIZE;

            // Replace with one offspring
            for (int j = 0; j < num_of_vertex; ++j)
                solutions[replaced_idx][j] = offsprings[i][j];

            // Recalculate value of replaced offsprings
            values[replaced_idx] = evaluate(graph_data, solutions[replaced_idx]);
        }

        // Update best, worst case solution
        min_avg_max = get_min_avg_max_from_vector(values, POPULATION_SIZE);

        best_solution_index = min_avg_max.max_idx;
        best_value = values[best_solution_index];

        double avg_of_values = min_avg_max.avg_value;
        printf("%.2f, %d, %d, %.2f\n", time_spent, generation, best_value, avg_of_values);
    }

    // Write output file
    write_out_file(out_file_name, solutions[best_solution_index], num_of_vertex);

    for(int i = 0; i < POPULATION_SIZE; i++)
        free(solutions[i]);
    SAFE_FREE(solutions);

    for(int i = 0; i < OFFSPRING_SIZE; i++)
        free(offsprings[i]);
    SAFE_FREE(offsprings);

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
