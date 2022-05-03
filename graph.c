
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