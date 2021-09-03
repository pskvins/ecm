//
// Created by 박석환 on 2021/09/02.
//

#include <algorithm>
#include <gsl/gsl_linalg.h>
#include <omp.h>
#include "newick.h"

void calculate_derivative(gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, newick_start *start, double *gradient) {
    gsl_matrix *dxepon_matrix[64 * 63 / 2];
    for (int num = 0; num < 64 * 63 / 2; num++) {
        dxepon_matrix[num] = gsl_matrix_alloc(64, 64);
    }
    #pragma omp parallel for schedule(dynamic)
    for (size_t num = 0; num < 64 * 63 / 2; num++) {

    }
}