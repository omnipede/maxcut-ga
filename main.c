#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "graph.h"
#include "util.h"

struct Graph read_in_file(char* filename);
void write_out_file(char* file_name, const int* vector, int vector_size);
void DFS(int vertex, int num_of_vertex, int* visited, struct Graph graph_data, int* mapping);
void BFS(int vertex, int num_of_vertex, int* visited, struct Graph graph_data, int* mapping);
void reorder(struct Graph graph_data, int** reordered_graph,int* mapping);

int newIdx=0, front=0, rear=0;

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
    double CROSSOVER_THRESHOLD = 0.249; // x10 (0 ~ 1)

    // Init randomness
    srand(time(NULL));

    // Read graph data
    struct Graph graph_data = read_in_file(in_file_name);
    //int** edges = graph_data.edges;
    int num_of_vertex = graph_data.num_of_vertex;

    // Reordering 
    int **reordered_graph = (int**) malloc(sizeof(int*) * num_of_vertex);
    for (int i = 0; i < num_of_vertex; ++i) {
        reordered_graph [i] = (int *) malloc(sizeof(int) * num_of_vertex);
        for (int j = 0; j < num_of_vertex; ++j)
            reordered_graph [i][j] = 0;
    }
    int* mapping = (int*)calloc(num_of_vertex,sizeof(int));
    reorder(graph_data, reordered_graph, mapping);
    
    struct Graph reordered_graph_data= init_graph(num_of_vertex, graph_data.num_of_edges, reordered_graph);
    int** edges = reordered_graph_data.edges;

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
        //values[i] = evaluate(graph_data, solutions[i]);
        values[i] = evaluate(reordered_graph_data, solutions[i]);
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

        // Selection
        int idx_of_mother = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);
        int idx_of_father = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);

        for (int i = 0; i < 10 && idx_of_mother == idx_of_father; ++i)
            idx_of_father = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);

        // Crossover
        int* child = (int*) malloc(sizeof (int) * num_of_vertex);
        for (int i = 0; i < num_of_vertex; ++i) {
            double r = rand() / (double)RAND_MAX;
            child[i] = r < CROSSOVER_THRESHOLD
                    ? solutions[idx_of_mother][i]
                    : solutions[idx_of_father][i];
        }

        // Mutation
        int maximum_mutation_count = (int)(num_of_vertex * 0.03);
        int mutated_count = rand() % maximum_mutation_count + 1;
        for (int i = 0; i < mutated_count; i++) {
            int mutated_index = rand() % num_of_vertex;
            child[mutated_index] = !child[mutated_index];
        }

        // Do local optimization
        //local_opt(graph_data, child);
        local_opt(reordered_graph_data, child);

        // Check equality
        int equal_solution_count = 0;
        for (int i = 0; i < POPULATION_SIZE; i++) {
            
            // Solution i 와 child 가 동일한지 확인한다.
            int is_equal = 1;
            for (int j = 0; j < num_of_vertex; j++) {
                if (solutions[i][j] != child[j]) {
                    is_equal = 0;
                    break;
                }
            }

            if (is_equal) {
                equal_solution_count += 1;
                break;
            }
        }

        if (equal_solution_count >= 1)
            continue;

        // Replace
        //int child_value = evaluate(graph_data, child);
        int child_value = evaluate(reordered_graph_data, child);
        int idx_to_replace = worst_solution_index;

        if (values[idx_of_father] >= child_value && values[idx_of_mother] < child_value)
            idx_to_replace = idx_of_mother;

        else if (values[idx_of_father] <= child_value && values[idx_of_mother] > child_value)
            idx_to_replace = idx_of_father;

        else if (values[idx_of_father] < child_value && values[idx_of_mother] < child_value) {

            int diff_with_mother = 0;
            int diff_with_father = 0;
            for (int i = 0; i < num_of_vertex; i++) {
                if (solutions[idx_of_mother][i] != child[i])
                    diff_with_mother += 1;
                if (solutions[idx_of_father][i] != child[i])
                    diff_with_father += 1;
            }

            // 두 부모 중에서 좀 더 유사한 부모와 변환
            idx_to_replace = (diff_with_mother > diff_with_father) 
                ? idx_of_father 
                : idx_of_mother;
        }

        // Replace with worst case
        for (int i = 0; i < num_of_vertex; ++i)
            solutions[idx_to_replace][i] = child[i];

        // Before update value, reduce sum of fitness
        sum_of_fitnesses -= fitnesses[idx_to_replace];

        // Update value of replaced solution
        values[idx_to_replace] = child_value;
        fitnesses[idx_to_replace] = (double)(values[idx_to_replace] - worst_value) + (best_value - worst_value) / (SELECTION_PRESSURE - 1.0);

        // Update sum of fitness
        sum_of_fitnesses += fitnesses[idx_to_replace];

        MACRO_FREE(child);

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

void DFS(int vertex, int num_of_vertex, int* visited, struct Graph graph_data, int* mapping){
    visited[vertex]=1;
    mapping[vertex]=newIdx;
    //printf("mapping: %d to %d\n",vertex,mapping[vertex]);
    newIdx++;

    for(int i=1;i<num_of_vertex;i++){
        if(graph_data.edges[vertex][i]!=0){
            if(visited[i]==0){
                DFS(i, num_of_vertex, visited, graph_data,mapping);
            }
        }
    }
}

void BFS(int vertex, int num_of_vertex, int* visited, struct Graph graph_data, int* mapping){
    int* queue = (int*)calloc(num_of_vertex,sizeof(int));
    visited[vertex]=1;
    queue[rear++]=vertex;

    while(1){
        if(rear<=front)break;
        
        vertex=queue[front++];
        visited[vertex]=1;
        mapping[vertex]=newIdx;
        printf("mapping: %d to %d\n",vertex,mapping[vertex]);
        newIdx++;

        for(int i=0;i<num_of_vertex;i++){
            if(graph_data.edges[vertex][i]!=0){
                if(visited[i]==0){
                    queue[rear++]=i;
                    visited[i]=1;
                }
            }
        }
    }
}

void reorder(struct Graph graph_data, int** reordered_graph,int* mapping){
    int num_of_vertex=graph_data.num_of_vertex;
    int* visited = (int*)calloc(num_of_vertex,sizeof(int));

    for(int i=0;i<num_of_vertex;i++){
        if(visited[i]==0){
            printf("start with vertex not visited: %d\n",i);
            DFS(i, graph_data.num_of_vertex, visited, graph_data, mapping);
        }
    }
    //BFS(0, graph_data.num_of_vertex, visited, graph_data, mapping);
    for(int i=0;i<num_of_vertex;i++){
        for(int j=0;j<num_of_vertex;j++){
            if(graph_data.edges[i][j]!=0){
                reordered_graph[mapping[i]][mapping[j]]=graph_data.edges[i][j];
            }
        }
    }
}