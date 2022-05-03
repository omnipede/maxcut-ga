
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
    // Variables
    int numOf0s = 0;
    int numOf1s = 0;
    int * s0 = NULL;
    int * s1 = NULL;
    int idxOfS0 = 0;
    int idxOfS1 = 0;
    int sum = 0;

    int sol_length = graph.num_of_vertex;
    int** weight_table = graph.weight_table;

    for(int i = 0; i < sol_length; i++) {
        if(sol[i] == 0)
            numOf0s += 1;
        else
            numOf1s += 1;
    }

    s0 = (int*)malloc(numOf0s*sizeof(int));
    s1 = (int*)malloc(numOf1s*sizeof(int));

    // read a gene of a solution in turn
    for(int i = 0; i < sol_length; i++) {
        if(sol[i] == 0) {
            s0[idxOfS0] = i;
            idxOfS0++;
        }

        if(sol[i] == 1) {
            s1[idxOfS1] = i;
            idxOfS1++;
        }
    }

    // Calculate the sum of weights by adding the weight of edges between s0 and s1
    for(int i = 0; i < numOf0s; i++) {
        for(int j = 0; j < numOf1s; j++) {
            sum += weight_table[s0[i]][s1[j]];
        }
    }

    // release memory
    free(s0);
    free(s1);

    return sum;
}

void clear_graph(struct Graph graph) {
    if (graph.weight_table) {
        for (int i = 0; i < graph.num_of_vertex; ++i) {
            if (graph.weight_table[i]) {
                free(graph.weight_table[i]);
                graph.weight_table[i] = NULL;
            }
        }
        free(graph.weight_table);
        graph.weight_table = NULL;
    }
}