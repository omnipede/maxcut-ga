#include <time.h>
#include <stdio.h>
#include <stdlib.h>
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

int newIdx=0;

void DFS(int vertex, int num_of_vertex, int* visited, struct Graph graph_data, int* mapping){
    visited[vertex]=1;
    mapping[vertex]=newIdx;
    newIdx++;

    for(int i=1;i<num_of_vertex;i++){
        if(graph_data.edges[vertex][i]!=0){
            if(visited[i]==0){
                DFS(i, num_of_vertex, visited, graph_data,mapping);
            }
        }
    }
}

void reorder(struct Graph graph_data, int** reordered_graph){
    int num_of_vertex=graph_data.num_of_vertex;
    int* visited = (int*)calloc(num_of_vertex,sizeof(int));
    int* mapping = (int*)calloc(num_of_vertex,sizeof(int));
    
    for(int i=0;i<num_of_vertex;i++){
        if(visited[i]==0){
            printf("start with vertex not visited: %d\n",i);
            DFS(i, graph_data.num_of_vertex, visited, graph_data, mapping);
        }
    }
    for(int i=0;i<num_of_vertex;i++){
        for(int j=0;j<num_of_vertex;j++){
            if(graph_data.edges[i][j]!=0&&reordered_graph[mapping[j]][mapping[i]]==0){
                reordered_graph[mapping[i]][mapping[j]]=graph_data.edges[i][j];
            }
        }
    }
}

int main(int argc, char *argv[]) {

    char* in_file_name = argv[1];
    char* out_file_name = argv[2];

    // Hyper parameters
    int POPULATION_SIZE = 500; // 100, 200, 300, 400, ..., 1000
    double SELECTION_PRESSURE = 3.4; // x10 (3 ~ 4)
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
    reorder(graph_data, reordered_graph);
    
    struct Graph reordered_graph_data= init_graph(num_of_vertex, graph_data.num_of_edges, reordered_graph);
    int** edges = reordered_graph_data.edges;

// DFS reordering 확인
/*
    int temp=0;
    for (int i = 0; i < num_of_vertex; ++i) {
        for (int j = 0; j < num_of_vertex; ++j){
            if(edges[i][j]!=0){
                temp++;
                printf("[%d][%d]:%d\n",i,j,edges[i][j]);
            }
        }
    }
    printf("%d: num of edge\n",temp);
    return 0;
*/
    

/*
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
        values[i] = evaluate(reordered_graph_data, solutions[i]);
    }
*/
    // 해집합 초기화시 중복해 제거 
    // Init population data
    int** solutions = (int**)malloc(POPULATION_SIZE * sizeof(int*));
    double* fitnesses = (double*)malloc(POPULATION_SIZE * sizeof(double));
    int* values = (int*)malloc(POPULATION_SIZE * sizeof(int));
    for(int i = 0; i < POPULATION_SIZE; ++i) {
        solutions[i] = (int*)malloc(num_of_vertex * sizeof(int));
        // Generate the set of solutions with random values
        while(1){
            int isUniq1=1, isUniq2=1;
            for(int j = 0; j < num_of_vertex; ++j)
                solutions[i][j] = rand() % 2;
            for(int p=0;p<i-1;p++){
                // 이전에 만들어진 해 중 중복해 있는지 검사 
                for(int j=0;j<num_of_vertex;++j){
                    if(solutions[i][j]!=solutions[p][j]){
                        isUniq1=1;
                        break;
                    } // 다르면 
                    else isUniq1=0; // 같으면 
                }
                // 반전시켜 중복해 있는지 검사
                for(int j=0;j<num_of_vertex;++j){
                    if((!solutions[i][j])!=solutions[p][j]){
                        isUniq2=1;
                        break;
                    }
                    else isUniq2=0;
                }
            }
            if(isUniq1==1&&isUniq2==1){
                fitnesses[i] = 0;
                values[i] = evaluate(reordered_graph_data, solutions[i]);
                break;
            }   
        }
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
    double EXECUTION_TIME = num_of_vertex/6;
    clock_t start = clock();
    while(1) {
        clock_t now = clock();
        double time_spent = (double)(now - start) / CLOCKS_PER_SEC;

        if (time_spent > EXECUTION_TIME)
            break;
/*
        // Selection
        int idx_of_mother = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);
        int idx_of_father = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);

        for (int i = 0; i < 10 && idx_of_mother == idx_of_father; ++i)
            idx_of_father = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);


        // Uniform Crossover
        int* child = (int*) malloc(sizeof (int) * num_of_vertex);
        for (int i = 0; i < num_of_vertex; ++i) {
            double r = rand() / (double)RAND_MAX;
            child[i] = r < CROSSOVER_THRESHOLD
                    ? solutions[idx_of_mother][i]
                    : solutions[idx_of_father][i];
        }
*/
/*
        // Cutpoint Crossover
        int CUT_NUM=5;
        int* child = (int*) malloc(sizeof (int) * num_of_vertex);
        int* cutpoint = (int*) calloc(CUT_NUM, sizeof (int));
        for(int i=0;i<CUT_NUM;i++){
            cutpoint[i]=rand()%num_of_vertex;
        }
        for(int i=0;i<CUT_NUM-1;i++){
            for(int j=1;j<CUT_NUM;j++){
                if(cutpoint[j]>cutpoint[i]){
                    int temp=cutpoint[j];
                    cutpoint[j]=cutpoint[i];
                    cutpoint[i]=temp;
                }
            }
        }
        
        for(int i=0;i<CUT_NUM-1;i++){
            double r = rand() / (double)RAND_MAX;
            if(r < CROSSOVER_THRESHOLD){
                for(int j=cutpoint[i];j<cutpoint[i+1];j++){
                    child[j] = solutions[idx_of_mother][j];
                }
            }else{
                for(int j=cutpoint[i];j<cutpoint[i+1];j++){
                    child[j] = solutions[idx_of_father][j];
                }
            }
        }
*/

        // Mutation
        /*
        int mutated_index = rand() % num_of_vertex;
        child[mutated_index] = !child[mutated_index];
        */

        // 대체 전에 중복해 있는지 검사, 중복해면 뮤테이션 
        int* child = (int*) malloc(sizeof (int) * num_of_vertex);
        while(1){
            // Selection
            int idx_of_mother = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);
            int idx_of_father = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);

            for (int i = 0; i < 10 && idx_of_mother == idx_of_father; ++i)
                idx_of_father = select_one_from_vector(fitnesses, POPULATION_SIZE, sum_of_fitnesses);


            // Uniform Crossover
            for (int i = 0; i < num_of_vertex; ++i) {
                double r = rand() / (double)RAND_MAX;
                child[i] = r < CROSSOVER_THRESHOLD
                        ? solutions[idx_of_mother][i]
                        : solutions[idx_of_father][i];
            }
            int isUniq1=1, isUniq2=1;
            int mutated_index = rand() % num_of_vertex;
            child[mutated_index] = !child[mutated_index];
            for(int p=0;p<POPULATION_SIZE;p++){
                // 이전에 만들어진 해 중 중복해 있는지 검사 
                for(int j=0;j<num_of_vertex;++j){
                    if(child[j]!=solutions[p][j]){
                        isUniq1=1;
                        break;
                    } // 다르면 
                    else isUniq1=0; // 같으면 
                }
                // 반전시켜 중복해 있는지 검사
                for(int j=0;j<num_of_vertex;++j){
                    if((!child[j])!=solutions[p][j]){
                        isUniq2=1;
                        break;
                    }// 다르면 
                    else isUniq2=0;// 같으면 
                }
            }
            if(isUniq1==1&&isUniq2==1)break;  
        }


        // Replace with worst case
        for (int i = 0; i < num_of_vertex; ++i)
            solutions[worst_solution_index][i] = child[i];

        // Before update value, reduce sum of fitness
        sum_of_fitnesses -= fitnesses[worst_solution_index];

        // Update value of replaced solution
        values[worst_solution_index] = evaluate(reordered_graph_data, child);
        fitnesses[worst_solution_index] = (double)(values[worst_solution_index] - worst_value) + (best_value - worst_value) / (SELECTION_PRESSURE - 1.0);

        // Update sum of fitness
        sum_of_fitnesses += fitnesses[worst_solution_index];

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
    // TODO: 리오더링한 버텍스 다시 변환해 결과파일 작성하는 작업 필요
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
