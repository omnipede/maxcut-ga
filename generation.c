#include "generation.h"

/* Private methods */
int select_parent_idx(struct Generation generation);
void replace(struct Generation* generation, struct Solution child);

double fitness(struct Generation generation, int idx);
void update_best_and_worst_and_avg_sols(struct Generation* generation);
void update_sum_of_fitness(struct Generation* generation);
////

/**
 * Generation 초기화 메소드
 */
void init_generation(struct Generation* generation, int population_size, struct Graph graph) {

    // Create population with random genes
    generation->population_size = population_size;
    generation->population = (struct Solution*) malloc(sizeof(struct Solution) * population_size);

    // Evaluate each solution
    generation->values = (int*) malloc(sizeof (int) * population_size);
    generation->graph = graph;

    int solution_size = graph.num_of_vertex;
    for (int i = 0; i < population_size; ++i) {
        init_solution(&generation->population[i], solution_size);
        generation->values[i] = evaluate(generation->graph, generation->population[i].genes);
    }

    // Update required values
    generation->sum_of_fitness = 0.0;
    update_best_and_worst_and_avg_sols(generation);
    update_sum_of_fitness(generation);

    // Set generation count
    generation->num_generation = 0;
}

void evolve(struct Generation* generation) {

    // Select
    int parent1_idx = select_parent_idx(*generation);
    int parent2_idx = select_parent_idx(*generation);
    for (int i = 0; i < 100 && parent1_idx == parent2_idx; ++i)
        parent2_idx = select_parent_idx(*generation);

    struct Solution* parent1, *parent2;
    parent1 = &generation->population[parent1_idx];
    parent2 = &generation->population[parent2_idx];

    // Crossover
    // TODO static 하게 바꾸기
    struct Solution* child = (struct Solution*) malloc(sizeof (struct Solution));
    crossover(*parent1, *parent2, child);

    // Mutate
    mutate(child);

    // Replace
    replace(generation, *child);

    // Clear child
    clear_solution(child);

    generation->num_generation += 1;
}

void print_generation(struct Generation generation) {
    printf("%d th generation best value: %d, avg value: %f, worst value:  %d\n", generation.num_generation,  generation.best_value, generation.avg_value, generation.worst_value);
}

struct Solution find_best_solution(struct Generation generation) {
    int best_index = generation.best_solution_idx;
    return generation.population[best_index];
}

void clear_generation(struct Generation* generation) {
    if (generation->population) {
        for (int i = 0; i < generation->population_size; ++i)
            clear_solution(&generation->population[i]);
        free(generation->population);
        generation->population = NULL;
    }

    if (generation->values) {
        free(generation->values);
        generation->values = NULL;
    }
}

int select_parent_idx(struct Generation generation) {
    double pt = (double) rand() / ((double)RAND_MAX + 1.0) * generation.sum_of_fitness;
    double sum = 0;

    for (int i = 0; i < generation.population_size; ++i) {
        sum += fitness(generation, i);
        if (pt < sum)
            return i;
    }

    return 0;
}

void replace(struct Generation* generation, struct Solution child) {
    // Replace with the worst solution
    int worst_idx = generation->worst_solution_idx;
    for (int i = 0; i < child.size; ++i)
        generation->population[worst_idx].genes[i] = child.genes[i];
    generation->values[worst_idx] = evaluate(generation->graph, child.genes);
    update_best_and_worst_and_avg_sols(generation);
    update_sum_of_fitness(generation);
}

void update_best_and_worst_and_avg_sols(struct Generation* generation) {

    generation->best_solution_idx = 0;
    generation->best_value = generation->values[generation->best_solution_idx];

    generation->worst_solution_idx = 0;
    generation->worst_value = generation->values[generation->worst_solution_idx];

    int sum = 0;
    for (int i = 0; i < generation->population_size; ++i) {
        if (generation->values[i] > generation->best_value) {
            generation->best_value = generation->values[i];
            generation->best_solution_idx = i;
        }

        if (generation->values[i] < generation->worst_value) {
            generation->worst_value = generation->values[i];
            generation->worst_solution_idx = i;
        }
        sum += generation->values[i];
    }

    generation->avg_value = sum / (double) generation->population_size;
}

void update_sum_of_fitness(struct Generation* generation) {

    double sum = 0;
    for (int i = 0; i < generation->population_size; ++i)
        sum += fitness(*generation, i);

    generation->sum_of_fitness = sum;
}

double fitness(struct Generation generation, int idx) {
    int ci = generation.values[idx];
    int cw = generation.worst_value;
    int cb = generation.best_value;

    return (double)(ci - cw) + (cb - cw) / (SELECTION_PRESSURE - 1.0);
}