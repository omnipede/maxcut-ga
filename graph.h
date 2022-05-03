

#ifndef MAXCUT_OPT_GRAPH2_H
#define MAXCUT_OPT_GRAPH2_H

#include <stdlib.h>
#include "util.h"

struct Graph {
    int num_of_vertex;
    int num_of_edges;
    int **weight_table;
};

struct Graph init_graph(int v, int e, int** table);
int evaluate(struct Graph, int* sol);

#endif //MAXCUT_OPT_GRAPH2_H
