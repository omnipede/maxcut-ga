
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

void local_opt(struct Graph graph, int* solution) {

    int num_of_vertex = graph.num_of_vertex;
    int* gains = (int*)malloc(sizeof(int) * num_of_vertex);
    int* locked = (int*)malloc(sizeof(int) * num_of_vertex);
    int* vis = (int*) malloc(sizeof(int) * num_of_vertex);
    int* gvis = (int*) malloc(sizeof(int) * num_of_vertex);

    while(1) {

        // Compute gains of all vertices
        for (int i = 0; i < num_of_vertex; i++) {
            int in_degree = 0;
            int out_degree = 0;

            int* weights = graph.edges[i];
            
            for (int j = 0; j < num_of_vertex; j += 1) {
                if (solution[i] == solution[j])
                    in_degree += weights[j];
                else
                    out_degree += weights[j];
            }

            gains[i] = in_degree - out_degree;
        }
        
        // Initialize set Q
        for (int i = 0; i < num_of_vertex; i += 1)
            locked[i] = 0;

        // Iterate for V - 1 times
        int vis_idx = 0;
        for (int n = 0; n < num_of_vertex - 1; n++) {
            
            // Find vertex who have maximal gain
            int vi = 0;
            int maximal_gain = INT_MIN;
            for (int i = 0; i < num_of_vertex; i++) {
                if (locked[i])
                    continue;
                
                if (gains[i] > maximal_gain) {
                    maximal_gain = gains[i];
                    vi = i;
                }
            }

            // Lock vertex vi
            locked[vi] = 1;

            // Insert vi into list
            vis[vis_idx] = vi;
            gvis[vis_idx] = maximal_gain;
            vis_idx += 1;

            for (int v = 0; v < num_of_vertex; v++) {
                
                // If v is in set Q, pass
                if (locked[v])
                    continue;

                // If v is in same set with vi
                if (solution[v] == solution[vi])
                    gains[v] -= 2 * graph.edges[vi][v];
                else
                    gains[v] += 2 * graph.edges[vi][v];
            }
        }

        // Choose k that maximize sum of gains from i = 0 to i = k
        int max_sum_of_gains = INT_MIN;
        int K = 0;

        for (int k = 0; k < num_of_vertex - 1; k++) {
            int sum_of_gains = 0;
            for (int i = 0; i <= k; i++)
                sum_of_gains += gvis[i];
            
            if (sum_of_gains > max_sum_of_gains) {
                max_sum_of_gains = sum_of_gains;
                K = k;
            }
        }

        if (max_sum_of_gains <= 0)
            break;

        for (int i = 0; i <= K; i++) {
            int vi = vis[i];
            solution[vi] = !solution[vi];
        }
    }

    MACRO_FREE(vis);
    MACRO_FREE(gvis);
    MACRO_FREE(locked);
    MACRO_FREE(gains);
}