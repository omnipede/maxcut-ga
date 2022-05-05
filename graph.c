
#include "graph.h"

/**
 * 그래프 데이터 초기화
 * @param v 정점 개수
 * @param e 간선 개수
 * @param table 간선 정보
 * @return
 */
struct Graph init_graph(int v, int e, int** table) {

    struct Graph graph_data = {
            v, e, table
    };

    return graph_data;
}

int evaluate(struct Graph graph, int* sol) {

    int solution_size = graph.num_of_vertex;

    // Split a graph vertexes into two sets
    int * first_set = (int*)malloc(solution_size*sizeof(int));
    int * second_set = (int*)malloc(solution_size*sizeof(int));

    // Put 0 in first_set and put 1 in second_set
    int num_of_0 = 0, num_of_1 = 0;
    for(int i = 0; i < solution_size; i++) {
        if(sol[i] == 0) {
            first_set[num_of_0] = i;
            num_of_0 += 1;
        }

        if(sol[i] == 1) {
            second_set[num_of_1] = i;
            num_of_1 += 1;
        }
    }

    // Calculate weight between two sets
    int sum = 0;
    int** edges = graph.edges;
    for(int i = 0; i < num_of_0; i++)
        for(int j = 0; j < num_of_1; j++) {
            int from = first_set[i];
            int to = second_set[j];
            sum += edges[from][to];
        }

    // release memory
    if (first_set)
        free(first_set);
    if (second_set)
        free(second_set);

    return sum;
}