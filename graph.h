//
// Created by 서현규 on 2022/05/03.
//

#ifndef MAXCUT_OPT_GRAPH2_H
#define MAXCUT_OPT_GRAPH2_H

struct Graph {
    int num_of_vertex;
    int num_of_edges;
    int **weight_table;
};

struct Graph init_graph(int v, int e, int** table);

#endif //MAXCUT_OPT_GRAPH2_H
