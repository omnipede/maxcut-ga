#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include "graph.h"

#define SAFE_FREE(a) if(a){free(a); a=NULL;} // for dynamic memory release

// Get the sum of weights of a solution
int get_sum_of_weight(int *sol, int sol_length, int **weight_table);

// Get the index of parent for Roulette Wheel selection from Solutions
int get_idx_of_parents(double *fitness, int sol_size, double fitnessSum);

struct Graph read_in_file(char* filename);

void write_out_file(char* filename);

int main(int argc, char *argv[])
{
    // Hyper parameters
    int POPULATION_SIZE = 400;
    double SELECTION_PRESSURE = 3.0;
    double CROSSOVER_THRESHOLD = 0.90;

    // Init randomness
    srand(time(NULL));

    // Read graph data
    struct Graph graph_data = read_in_file("../maxcut.in");
    int** weight_table = graph_data.weight_table;
    int num_of_vertex = graph_data.num_of_vertex;
    int num_of_edges = graph_data.num_of_edges;

    // Init population data
    int** solutions = (int**)malloc(POPULATION_SIZE * sizeof(int*));
    double* fitnesses = (double*)malloc(POPULATION_SIZE * sizeof(double));
    int* values = (int*)malloc(POPULATION_SIZE * sizeof(int));
    for(int i = 0; i < POPULATION_SIZE; i++) {
        solutions[i] = (int*)malloc(num_of_vertex * sizeof(int));
        // Generate the set of solutions with random values
        for(int j = 0; j < num_of_vertex; j++)
            solutions[i][j] = (int)(((double)rand() / ((double)(RAND_MAX) + 1.0)) * 2);
        fitnesses[i] = 0;
        values[i] = get_sum_of_weight(solutions[i], num_of_vertex, weight_table);
    }

    // find worst and best values
    int worst_solution_index = 0, worst_value = values[worst_solution_index];
    int best_solution_index = 0, best_value = values[best_solution_index];
    double sum_of_fitnesses = 0;

    for (int i = 0; i < POPULATION_SIZE; ++i) {
        if (values[i] > best_value) {
            best_value = values[i];
            best_solution_index = i;
        }

        if (values[i] < worst_value) {
            worst_value = values[i];
            worst_solution_index = i;
        }
    }

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

        if (time_spent > 10)
            break;

        // Selection
        int idx_of_mother = get_idx_of_parents(fitnesses, POPULATION_SIZE, sum_of_fitnesses);
        int idx_of_father = get_idx_of_parents(fitnesses, POPULATION_SIZE, sum_of_fitnesses);

        while(idx_of_mother == idx_of_father)
            idx_of_father = get_idx_of_parents(fitnesses, POPULATION_SIZE, sum_of_fitnesses);

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
        child[mutated_index] = child[mutated_index] == 1 ? 0 : 1;

        // Replace with worst case
        for (int i = 0; i < num_of_vertex; ++i)
            solutions[worst_solution_index][i] = child[i];

        // Update value of replaced solution
        values[worst_solution_index] = get_sum_of_weight(child, num_of_vertex, weight_table);

        SAFE_FREE(child);

        // Update best, worst case solution
        worst_solution_index = 0;
        worst_value = values[worst_solution_index];
        best_solution_index = 0;
        best_value = values[best_solution_index];

        int sum_of_values = 0;
        for (int i = 0; i < POPULATION_SIZE; ++i) {
            if (values[i] > best_value) {
                best_value = values[i];
                best_solution_index = i;
            }

            if (values[i] < worst_value) {
                worst_value = values[i];
                worst_solution_index = i;
            }

            sum_of_values += values[i];
        }

        double avg_of_values = sum_of_values / (double) POPULATION_SIZE;

        // Update sum of fitness
        sum_of_fitnesses = 0;
        for (int i = 0; i < POPULATION_SIZE; ++i) {
            int ci = values[i];
            int cw = worst_value;
            int cb = best_value;

            fitnesses[i] = (double)(ci - cw) + (cb - cw) / (SELECTION_PRESSURE - 1.0);
            sum_of_fitnesses += fitnesses[i];
        }

        printf("%.2f, %d, %.2f\n", time_spent, best_value, avg_of_values);
    }

    // Write output file
    FILE* out_file = fopen("../maxcut.out", "w");
    int flag = solutions[best_solution_index][0];
    for (int i = 0; i < num_of_vertex; ++i)
        // Warning: start index from '1'
        if (solutions[best_solution_index][i] == flag)
            fprintf(out_file, "%d ", i + 1);
    fclose(out_file);

    for(int i = 0; i < POPULATION_SIZE; i++)
        free(solutions[i]);
    SAFE_FREE(solutions);

    for(int i = 0; i < num_of_vertex; i++)
        free(weight_table[i]);
    SAFE_FREE(weight_table);

    return 0;
}


// Get the sum of weights of a solution
// input: parents[i], length of parents[i], weight_table
// output: sum of weights of parents[i]
int get_sum_of_weight(int * sol, int sol_length, int ** weight_table)
{
    // Variables
    int numOf0s = 0;
    int numOf1s = 0;
    int * s0 = NULL;
    int * s1 = NULL;
    int idxOfS0 = 0;
    int idxOfS1 = 0;
    int sum = 0;

    // Find # of 0s in a solution
    //printf("\nSOL: ");
    for(int i = 0; i < sol_length; i++)
    {
        if(sol[i] == 0)
            numOf0s++;
        //printf("%d", sol[i]);
    }
    numOf1s = sol_length - numOf0s;
    //printf("\n0s: %d, 1s: %d\n", numOf0s, numOf1s);

    // split a solution into two groups, s0 & s1
    // S0 for the index of 0s in a solution
    // S1 for the index of 1s in a solution
    s0 = (int*)malloc(numOf0s*sizeof(int));
    s1 = (int*)malloc(numOf1s*sizeof(int));
    // read a gene of a solution in turn
    // if 0, then write it to s0
    // else, then write it to s1
    for(int i = 0; i < sol_length; i++)
    {
        if(sol[i] == 0)
        {
            s0[idxOfS0] = i;
            //printf("V[%d]=%d->s0[%d]=%d\n", i, sol[i], idxOfS0, s0[idxOfS0]);
            idxOfS0++;
        }
    }
    for(int i = 0; i < sol_length; i++)
    {
        if(sol[i] == 1)
        {
            s1[idxOfS1] = i;
            //printf("V[%d]=%d->s1[%d]=%d\n", i, sol[i], idxOfS1, s1[idxOfS1]);
            idxOfS1++;
        }
    }
    //printf("\n");

    // Calculate the sum of weights by adding the weight of edges between s0 and s1
    for(int i = 0; i < numOf0s; i++)
    {
        for(int j = 0; j < numOf1s; j++)
        {
            sum += weight_table[s0[i]][s1[j]];
            //printf("sum=%d\n", sum);
        }
    }

    // release memory
    SAFE_FREE(s0);
    SAFE_FREE(s1);

    return sum;
}

// Get the index of parent for Roulette Wheel selection from Solutions
// input: fitnesses, # of solutions, sum of fitnesses of solutions
// output: index of selected solution
int get_idx_of_parents(double *fitness, int sol_size, double fitnessSum)
{
    double point = (double)rand() / ((double)RAND_MAX + 1.0) * fitnessSum;
    //printf("fitness sum: %f, point: %.2f\n", fitnessSum, point);
    double sum = 0.0;

    for(int i = 0; i < sol_size; i++)
    {
        sum += fitness[i];
        if(point < sum)
        {
            return i;
        }
    }

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
