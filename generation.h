#ifndef MAXCUT_OPT_GENERATION_H
#define MAXCUT_OPT_GENERATION_H

#include "solution.h"
#include "graph.h"

#define SELECTION_PRESSURE 3.0

struct Generation {

    int population_size;
    struct Solution* population;

    struct Graph graph;
    int* values;

    int best_value;
    int best_solution_idx;
    int worst_value;
    int worst_solution_idx;
    double avg_value;

    double sum_of_fitness;
    int num_generation;
};

void init_generation(struct Generation* generation, int population_size, struct Graph graph);
void evolve(struct Generation* generation);
void print_generation(struct Generation generation);
struct Solution find_best_solution(struct Generation generation);
void clear_generation(struct Generation* generation);

#endif //MAXCUT_OPT_GENERATION_H
