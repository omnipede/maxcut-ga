
#ifndef MAXCUT_OPT_SOLUTION_H
#define MAXCUT_OPT_SOLUTION_H

#include <stdlib.h>
#include <stdio.h>

#define CROSSOVER_THRESHOLD 0.90

struct Solution {
    int size;
    int* genes;
};

void init_solution(struct Solution* solution, int size);
void init_solution_with_genes(struct Solution* solution, int size, int* genes);

void crossover(struct Solution parent1, struct Solution parent2, struct Solution* child);
void mutate(struct Solution* sol);

void print_solution(struct Solution sol);
void clear_solution(struct Solution* solution);

#endif //MAXCUT_OPT_SOLUTION_H
