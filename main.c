#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "graph.h"
#include "solution.h"
#include "generation.h"

struct Graph read_in_file(char* filename);
void write_out_file(char* filename, struct Generation generation);

int main(int argc, char *argv[]) {
    // Init randomness
    srand(time(NULL) + (unsigned)getpid());

    // Read input file
    struct Graph graph = read_in_file("../maxcut.in");
    struct Generation generation;
    init_generation(&generation, 400, graph);

    clock_t start = clock();
    while (1) {
        clock_t now = clock();
        double time_spent = (double)(now - start) / CLOCKS_PER_SEC;

        if (time_spent > 10)
            break;

        evolve(&generation);
        print_generation(generation);
    }

    // Write output file
    write_out_file("../maxcut.out", generation);

    clear_generation(&generation);
    clear_graph(graph);

    return 0;
}

struct Graph read_in_file(char* filename) {
    int num_of_vertex, num_of_edges;
    int **graph = NULL;
    FILE* in_file = fopen(filename, "r");
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

void write_out_file(char* filename, struct Generation generation) {
    // Write output file
    FILE* out_file = fopen(filename, "w");
    struct Solution sol = find_best_solution(generation);
    int flag = sol.genes[0];
    for (int i = 0; i < sol.size; ++i)
        if (sol.genes[i] == flag)
            fprintf(out_file, "%d ", i + 1);
    fclose(out_file);
}
