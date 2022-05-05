
#ifndef MAXCUT_OPT_UTIL_H
#define MAXCUT_OPT_UTIL_H

#define SAFE_FREE(a) if(a){free(a); a=NULL;}

#include <stdlib.h>
#include <stdio.h>

int select_one_from_vector(const double* vector, int vector_size, double sum_of_vector);
void select_two_from_vector(const int* vector, int vector_size, int* buffer);

struct MinAvgMax {
    int min_idx;
    double avg_value;
    int max_idx;
};

struct MinAvgMax get_min_avg_max_from_vector(const int* vector, int vector_size);

#endif //MAXCUT_OPT_UTIL_H
