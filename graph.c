
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

    // Find # of 0s in a solution
    int numOf0s = 0, numOf1s = 0;
    int sol_length = graph.num_of_vertex;
    for(int i = 0; i < sol_length; i++){
        if(sol[i] == 0)
            numOf0s += 1;
        else
            numOf1s += 1;
    }

    // split a solution into two groups, s0 & s1
    int * s0 = (int*)malloc(numOf0s*sizeof(int));
    int * s1 = (int*)malloc(numOf1s*sizeof(int));

    // if 0, then write it to s0
    // else, then write it to s1
    int idx_of_s0 = 0;
    int idx_of_s1 = 0;
    for(int i = 0; i < sol_length; i++) {
        if(sol[i] == 0) {
            s0[idx_of_s0] = i;
            idx_of_s0++;
        }

        if(sol[i] == 1) {
            s1[idx_of_s1] = i;
            idx_of_s1++;
        }
    }

    int sum = 0;
    int** weight_table = graph.weight_table;
    for(int i = 0; i < numOf0s; i++) {
        for(int j = 0; j < numOf1s; j++) {
            sum += weight_table[s0[i]][s1[j]];
        }
    }

    // release memory
    SAFE_FREE(s0);
    SAFE_FREE(s1);

    return sum;
}