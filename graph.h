

#ifndef MAXCUT_OPT_GRAPH2_H
#define MAXCUT_OPT_GRAPH2_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"

struct Graph {
    int num_of_vertex;
    int num_of_edges;
    int **edges;
};

struct Graph init_graph(int v, int e, int** table);
int evaluate(struct Graph, int* sol);
void local_opt(struct Graph graph, int* solution);

#endif //MAXCUT_OPT_GRAPH2_H
