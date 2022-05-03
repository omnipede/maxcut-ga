#include "solution.h"

/**
 * Solution (해) 생성자
 * @param solution
 * @param size
 */
void init_solution(struct Solution* solution, int size) {
    int* genes = (int*) malloc (sizeof(int) * size);
    for (int i = 0; i < size; ++i)
        genes[i] = (int)(((double) rand() / ((double)(RAND_MAX) + 1.0)) * 2);

    solution->genes = genes;
    solution->size = size;
}

void init_solution_with_genes(struct Solution* solution, int size, int* genes) {
    solution->genes = genes;
    solution->size = size;
}

void crossover(struct Solution parent1, struct Solution parent2, struct Solution* child) {
    if (parent1.size != parent2.size)
        return;

    int* child_genes = (int*) malloc (sizeof (int) * parent1.size);
    for (int i = 0; i < parent1.size; ++i) {
        double r = rand() / ((double)RAND_MAX + 1.0);
        child_genes[i] = r < CROSSOVER_THRESHOLD
                   ? parent1.genes[i]
                   : parent2.genes[i];
    }

    child->size = parent1.size;
    child->genes = child_genes;

//    init_solution_with_genes(child, parent1.size, child_genes);
}

void mutate(struct Solution* sol) {
    int mutated_index = (int) (rand() / ((double)RAND_MAX + 1.0) * sol->size) % sol->size;
    sol->genes[mutated_index] = sol->genes[mutated_index] == 1 ? 0 : 1;
}

void print_solution(struct Solution sol) {

    for (int i = 0; i < sol.size; ++i)
        printf("%d", sol.genes[i]);
    printf("\n");
}

void clear_solution(struct Solution* solution) {
    if (solution) {
        if (solution->genes) {
            free(solution->genes);
            solution->genes = NULL;
        }
    }
}
