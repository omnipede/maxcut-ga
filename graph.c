
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

int isDuplicated(int* child, int** solutions, int POPULATION_SIZE, int num_of_vertex){
    int isUniq1=1, isUniq2=1;
    for(int p=0;p<POPULATION_SIZE;++p){
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
    if(isUniq1==1&&isUniq2==1)return 0;
    else return 1;
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

void KL(struct Graph graph, int* solution) {

    int step = 0;
    while(1) {
        int num_of_vertex = graph.num_of_vertex;
        int* group_a = (int*) malloc( sizeof(int) * num_of_vertex);
        int* group_b = (int*) malloc( sizeof(int) * num_of_vertex);
        int group_a_idx = 0, group_b_idx = 0;

        for (int i = 0; i < num_of_vertex; i++) {
            if (solution[i] == 0)
                group_a[group_a_idx++] = i;
            else
                group_b[group_b_idx++] = i;
        }

        // Calc d_values for all vertices
        int* d_values = (int*) malloc(sizeof(int) * num_of_vertex);
        for (int i = 0; i < num_of_vertex; i++) {
            int in_degree = 0;
            int out_degree = 0;

            for (int j = 0; j < num_of_vertex; j += 1) {
                int weight = graph.edges[i][j];
                if (solution[i] == solution[j])
                    in_degree += weight;
                else
                    out_degree += weight;
            }

            d_values[i] = -1 * (out_degree - in_degree);
        }

        int* visited = (int*) malloc(sizeof(int) * num_of_vertex);
        for (int i = 0; i < num_of_vertex; i++)
            visited[i] = 0;

        int* gv = (int*) malloc(sizeof(int) * num_of_vertex);
        int* ga = (int*) malloc(sizeof(int) * num_of_vertex);
        int* gb = (int*) malloc(sizeof(int) * num_of_vertex);
        int g_idx = 0;

        for (int n = 0; n < num_of_vertex / 2; n++) {

            // Find maximal gain
            int max_gain = -1;
            int target_a = 0, target_b = 0;

            for (int i = 0; i < group_a_idx; i ++ ) {
                int a = group_a[i];
                if (visited[a])
                    continue;

                for (int j = 0; j < group_b_idx; j++) {

                    int b = group_b[j];
                    if (visited[b])
                        continue;

                    int c_ab = graph.edges[a][b];
                    int gain = d_values[a] + d_values[b] - 2 * c_ab;

                    if (gain > max_gain) {
                        max_gain = gain;
                        target_a = a;
                        target_b = b;
                    } 
                }
            }

            // Update results
            gv[g_idx] = max_gain;
            ga[g_idx] = target_a;
            gb[g_idx] = target_b;
            g_idx += 1;

            // Mark as visited
            visited[target_a] = 1;
            visited[target_b] = 1;

            for (int i = 0; i < group_a_idx; i++) {
                int x = group_a[i];
                if (visited[x])
                    continue;

                int c_xa = graph.edges[x][target_a];
                int c_xb = graph.edges[x][target_b];
                d_values[x] += -1 * ( 2 * c_xa - 2 * c_xb );
            }

            for (int i = 0; i < group_b_idx; i++) {
                int y = group_b[i];
                if (visited[y])
                    continue;
                
                int c_yb = graph.edges[y][target_b];
                int c_ya = graph.edges[y][target_a];
                d_values[y] += -1 * ( 2 * c_yb - 2 * c_ya );
            }
        }

        int max_sum_of_gain = -1;
        int k_max = 0;

        for (int k = 1; k < g_idx - 1; k++) {
            int sum_of_gain = 0;
            for (int i = 0; i <= k; i++)
                sum_of_gain += gv[i];

            if (sum_of_gain > max_sum_of_gain) {
                k_max = k;
                max_sum_of_gain = sum_of_gain;
            }
        }

        if (max_sum_of_gain > 0) {
            for (int i = 0; i <= k_max; i++) {
                int a = ga[i];
                int b = gb[i];

                solution[a] = !solution[a];
                solution[b] = !solution[b];
            }
        }

        // free
        MACRO_FREE(gv);
        MACRO_FREE(ga);
        MACRO_FREE(gb);

        MACRO_FREE(visited);
        MACRO_FREE(d_values);

        MACRO_FREE(group_a);
        MACRO_FREE(group_b);

        if (max_sum_of_gain <= 0)
            break;

        if (step >= 1)
            break;
        step += 1;
    }
}

void local_opt(struct Graph graph, int* solution) {
    
    int before_K = 0;
    while(1) {

        int num_of_vertex = graph.num_of_vertex;

        // Compute gains of all vertices
        int* gains = (int*)malloc(sizeof(int) * num_of_vertex);
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
        int* locked = (int*)malloc(sizeof(int) * num_of_vertex);
        for (int i = 0; i < num_of_vertex; i += 1)
            locked[i] = 0;

        // Iterate for V - 1 times
        int* vis = (int*) malloc(sizeof(int) * num_of_vertex);
        int* gvis = (int*) malloc(sizeof(int) * num_of_vertex);
        int vis_idx = 0;
        for (int n = 0; n < num_of_vertex - 1; n++) {
            
            // Find vertex who have maximal gain
            int vi = 0;
            int maximal_gain = INT32_MIN;
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
        int max_sum_of_gains = INT32_MIN;
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

        if (max_sum_of_gains == 0 || K == before_K) {
            MACRO_FREE(vis);
            MACRO_FREE(gvis);
            MACRO_FREE(locked);
            MACRO_FREE(gains);
            break;
        }

        before_K = K;

        for (int i = 0; i <= K; i++) {
            int vi = vis[i];
            solution[vi] = !solution[vi];
        }

        MACRO_FREE(vis);
        MACRO_FREE(gvis);
        MACRO_FREE(locked);
        MACRO_FREE(gains);
    }
}