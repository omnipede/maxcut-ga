#include "util.h"

/**
 * Get one from vector on certain method. (Roulette wheel, etc ...)
 * @param vector
 * @param vector_size
 * @param sum_of_vector
 * @return index
 */
int select_one_from_vector(const double* vector, int vector_size, double sum_of_vector) {

    double sum = 0;
    double pt = (double) rand() / ((double) RAND_MAX + 1.0) * sum_of_vector;
    for (int i = 0; i < vector_size; ++i) {
        sum += vector[i];
        if (pt < sum)
            return i;
    }

    return 0;
}

void select_two_from_vector(const int* vector, int vector_size, int* buffer) {

    for (int i = 0; i < 50; ++i) {
        int temp = rand() % vector_size;

        if (vector[temp] >= vector[buffer[0]]) {
            buffer[0] = temp;
        }

        if (vector[temp] <= vector[buffer[0]]) {
            buffer[1] = temp;
        }
    }
}

/**
 * Get minimum index, average value, maximum index from vector
 * @param vector
 * @param vector_size
 * @return
 */
struct MinAvgMax get_min_avg_max_from_vector(const int* vector, int vector_size) {
    int worst_index = 0, best_index = 0;
    int worst_value = vector[worst_index], best_value = vector[best_index];

    int sum_of_values = 0;
    for (int i = 0; i < vector_size; ++i) {
        int v = vector[i];

        if (v > best_value) {
            best_value = v;
            best_index = i;
        }

        if (v < worst_value) {
            worst_value = v;
            worst_index = i;
        }

        sum_of_values += v;
    }

    double avg_value = sum_of_values / (double) vector_size;

    struct MinAvgMax ret = {
            worst_index, avg_value, best_index
    };

    return ret;
}