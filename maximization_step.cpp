//
// Created by 박석환 on 2021/09/02.
//

#include <algorithm>
#include <gsl/gsl_linalg.h>
#include <omp.h>
#include "newick.h"
#include "maximization_step.h"

void calculate_derivative(gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, newick_start *start, double *gradient) {
    gsl_matrix *dxepon_matrix[64 * 63 / 2];
    for (int num = 0; num < 64 * 63 / 2; num++) {
        dxepon_matrix[num] = gsl_matrix_alloc(64, 64);
    }
    double diagonal[64] = {0};
    double normalize;
    for (int row = 0; row < 64; row++) {
        for (int col = 0; col < 64; col++) {
            if (row != col) {
                diagonal[row] -= gsl_matrix_get(qmatrix, row, col) * codon_freq[col];
            }
        }
        normalize -= diagonal[row] * codon_freq[row];
    }
    #pragma omp parallel for schedule(dynamic)
    for (size_t diff = 0; diff < 64 * 63 / 2; diff++) {
        thread_local std::vector<newick_graph*> next_iteration = start->next;
        thread_local gsl_matrix *dx_matrix = gsl_matrix_alloc(64, 64);//freed
        thread_local gsl_matrix *F = gsl_matrix_alloc(64, 64);//freed
        thread_local size_t row_tar = num_to_coordinate.find(diff)->second.begin()->first;
        thread_local size_t col_tar = num_to_coordinate.find(diff)->second.begin()->second;
        gsl_matrix_set_all(dx_matrix, 0.0);
        while (next_iteration.size() != 1) {
            size_t size = next_iteration.size();
            for (size_t num = 0; num < size; num++) {
                for (size_t row = 0; row < 64; row++) {
                    for (size_t col = 0; col < 64; col++) {
                        if (row == col) {
                            gsl_matrix_set(dx_matrix, row, col,
                                           (2 * codon_freq[row_tar] * codon_freq[col_tar] * diagonal[row] - codon_freq[col_tar] * normalize) / pow(normalize, 2));
                        } else if ((row == row_tar && col == col_tar) || (row == col_tar && col == row_tar)) {
                            gsl_matrix_set(dx_matrix, row, col,
                                           (codon_freq[col] * normalize - 2 * codon_freq[row] * pow(codon_freq[col], 2) * gsl_matrix_get(qmatrix, row_tar, col_tar) / pow(normalize, 2)));
                        } else {
                            gsl_matrix_set(dx_matrix, row, col,
                                           -gsl_matrix_get(qmatrix, row, col) * 2 * codon_freq[row_tar] * codon_freq[col_tar] / pow(normalize, 2));
                        }
                    }
                }
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigen_inverse, dx_matrix, 0.0, dx_matrix);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dx_matrix, eigenvector, 0.0, dx_matrix);
                for (size_t row = 0; row < 64; row++) {
                    for (size_t col = 0; col < 64; col++) {
                        if (abs(gsl_vector_get(eigenvalue, row) - gsl_vector_get(eigenvalue, col)) < 1.0e-14) {
                            gsl_matrix_set(F, row, col,
                                           next_iteration[num]->branch_length * pow(M_E, gsl_vector_get(eigenvalue, row) * next_iteration[num]->branch_length));
                        } else {
                            gsl_matrix_set(F, row, col,
                                           (pow(M_E, gsl_vector_get(eigenvalue, row) * next_iteration[num]->branch_length)
                                           - pow(M_E, gsl_vector_get(eigenvalue, col) * next_iteration[num]->branch_length))
                                            / (gsl_vector_get(eigenvalue, row) - gsl_vector_get(eigenvalue, col)));
                        }
                    }
                }
                gsl_matrix_mul_elements(dx_matrix, F);
                gsl_matrix_free(F);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigenvector, dx_matrix, 0.0, dx_matrix);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dx_matrix, eigen_inverse, 0.0, dxepon_matrix[num]);
                gsl_matrix_free(dx_matrix);
            }
        }
    }
}