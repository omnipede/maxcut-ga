//
// Created by 서현규 on 2022/05/03.
//

#ifndef MAXCUT_OPT_GRAPH2_H
#define MAXCUT_OPT_GRAPH2_H

#include <stdlib.h>

struct Graph {
    int num_of_vertex;
    int num_of_edges;
    int **weight_table;
};

struct Graph init_graph(int v, int e, int** table);
int evaluate(struct Graph graph, int* sol);
void clear_graph(struct Graph graph);

#endif //MAXCUT_OPT_GRAPH2_H
