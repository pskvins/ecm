//
// Created by 박석환 on 2021/10/18.
//

#include <algorithm>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "process_maf.h"
#include <assert.h>
#include <omp.h>

void get_inverse_of_eigenvec_com(gsl_matrix_complex *eigenvector_com, gsl_matrix_complex *eigenvec_inverse_com) {
    gsl_matrix_complex_memcpy(eigenvec_inverse_com, eigenvector_com);
    gsl_permutation *p = gsl_permutation_alloc(64);//freed
    int here = 0;
    int *signum = &here;
    gsl_linalg_complex_LU_decomp(eigenvec_inverse_com, p, signum);
    gsl_linalg_complex_LU_invx(eigenvec_inverse_com, p);
    gsl_permutation_free(p);
}

gsl_matrix* calculate_expon_matrix_com(gsl_matrix_complex *eigenvector_com, gsl_matrix_complex *eigenvec_inverse_com, gsl_vector_complex *eigenvalue_com,
                                       gsl_matrix *eigenvector_temp, gsl_vector *eigenvalue_temp, double branch_length, gsl_vector *ck, gsl_vector *sk) {
    double element;
    for (int i = 0; i < 64; i++) {
        gsl_vector_set(eigenvalue_temp, i, pow(M_E, gsl_vector_complex_get(eigenvalue_com, i).dat[0] * branch_length));
        gsl_vector_set(ck, i, gsl_vector_get(eigenvalue_temp, i) * cos(gsl_vector_complex_get(eigenvalue_com, i).dat[1] * branch_length));
        gsl_vector_set(sk, i, gsl_vector_get(eigenvalue_temp, i) * sin(gsl_vector_complex_get(eigenvalue_com, i).dat[1] * branch_length));
    }
    for (int row = 0; row < 64; row++) {
        for (int col = 0; col < 64; col++) {
            element = 0.0;
            for (int i = 0; i < 64; i++) {
                if (gsl_vector_complex_get(eigenvalue_com, i).dat[1] == 0.0) {
                    element += gsl_matrix_complex_get(eigenvector_com, row, i).dat[0] *
                            gsl_vector_get(eigenvalue_temp, i) * gsl_matrix_complex_get(eigenvec_inverse_com, i, col).dat[0];
                } else if (gsl_vector_complex_get(eigenvalue_com, i).dat[1] > 0.0) {
                    element += 2 * (gsl_vector_get(ck, i) *
                            (gsl_matrix_complex_get(eigenvector_com, row, i).dat[0] * gsl_matrix_complex_get(eigenvec_inverse_com, i, col).dat[0] -
                            gsl_matrix_complex_get(eigenvector_com, row, i).dat[1] * gsl_matrix_complex_get(eigenvec_inverse_com, i, col).dat[1]) -
                            gsl_vector_get(sk, i) *
                            (gsl_matrix_complex_get(eigenvector_com, row, i).dat[0] * gsl_matrix_complex_get(eigenvec_inverse_com, i, col).dat[1] +
                            gsl_matrix_complex_get(eigenvector_com, row, i).dat[1] * gsl_matrix_complex_get(eigenvec_inverse_com, i, col).dat[0]));
                }
            }
            gsl_matrix_set(eigenvector_temp, row, col, element);
            if (element > 1 || element < 0) {
                std::cout << "Probability out of range : " << element << std::endl;
                exit(-1);
            }
        }
    }
    return eigenvector_temp;
}

void set_matrices_complex(newick_start *start, gsl_matrix_complex *eigenvector_com, gsl_matrix_complex *eigen_inverse_com,
                          gsl_vector_complex *eigenvalue_com, int *newick_order_max) {
    std::vector<newick_graph*> next_iterator = start->next;
    gsl_vector *eigenvalue_temp = gsl_vector_alloc(64);//freed
    gsl_matrix *eigenvector_temp = gsl_matrix_alloc(64, 64);//freed
    gsl_vector *ck = gsl_vector_alloc(64);
    gsl_vector *sk = gsl_vector_alloc(64);
    gsl_complex element;
    int newick_order = *newick_order_max;
    while (newick_order >= 0) {
        size_t size = next_iterator.size();
        for (size_t num = 0; num < size; num++) {
            if (next_iterator[num]->order == newick_order) {
                double branch_length = next_iterator[num]->branch_length;
                gsl_matrix_memcpy(next_iterator[num]->expon_matrix, calculate_expon_matrix_com(eigenvector_com, eigen_inverse_com, eigenvalue_com, eigenvector_temp, eigenvalue_temp, branch_length, ck, sk));
                for (int row = 0; row < 64; row++) {
                    for (int col = 0; col < 64; col++) {
                        if (row == col) {
                            element = gsl_complex_rect(branch_length * gsl_vector_get(ck, row), branch_length * gsl_vector_get(sk, row));
                            gsl_matrix_complex_set(next_iterator[num]->expon_eigen_com, row, col, element);
                        }
                        if (row > col) {
                            if ((abs(gsl_vector_complex_get(eigenvalue_com, row).dat[0] - gsl_vector_complex_get(eigenvalue_com, col).dat[0]) < 1.0e-13) &&
                            (abs(gsl_vector_complex_get(eigenvalue_com, row).dat[1] - gsl_vector_complex_get(eigenvalue_com, col).dat[1]) < 1.0e-13)) {
                                element = gsl_complex_rect(branch_length * gsl_vector_get(ck, row), branch_length * gsl_vector_get(sk, row));
                                gsl_matrix_complex_set(next_iterator[num]->expon_eigen_com, row, col, element);
                                gsl_matrix_complex_set(next_iterator[num]->expon_eigen_com, col, row, element);
                            } else {
                                element.dat[0] = ((gsl_vector_get(ck, row) - gsl_vector_get(ck, col)) * (gsl_vector_complex_get(eigenvalue_com, row).dat[0] - gsl_vector_complex_get(eigenvalue_com, col).dat[0]) +
                                                (gsl_vector_get(sk, row) - gsl_vector_get(sk, col)) * (gsl_vector_complex_get(eigenvalue_com, row).dat[1] - gsl_vector_complex_get(eigenvalue_com, col).dat[1])) /
                                                (pow(gsl_vector_complex_get(eigenvalue_com, row).dat[0] - gsl_vector_complex_get(eigenvalue_com, col).dat[0], 2) +
                                                pow(gsl_vector_complex_get(eigenvalue_com, row).dat[1] - gsl_vector_complex_get(eigenvalue_com, col).dat[1], 2));
                                element.dat[1] = ((gsl_vector_get(sk, row) - gsl_vector_get(sk, col)) * (gsl_vector_complex_get(eigenvalue_com, row).dat[0] - gsl_vector_complex_get(eigenvalue_com, col).dat[0]) -
                                                   (gsl_vector_get(ck, row) - gsl_vector_get(ck, col)) * (gsl_vector_complex_get(eigenvalue_com, row).dat[1] - gsl_vector_complex_get(eigenvalue_com, col).dat[1])) /
                                                  (pow(gsl_vector_complex_get(eigenvalue_com, row).dat[0] - gsl_vector_complex_get(eigenvalue_com, col).dat[0], 2) +
                                                   pow(gsl_vector_complex_get(eigenvalue_com, row).dat[1] - gsl_vector_complex_get(eigenvalue_com, col).dat[1], 2));
                                gsl_matrix_complex_set(next_iterator[num]->expon_eigen_com, row, col, element);
                                gsl_matrix_complex_set(next_iterator[num]->expon_eigen_com, col, row, element);
                            }
                        }
                    }
                }
                if (std::find(next_iterator.begin(), next_iterator.end(), next_iterator[num]->next) == next_iterator.end()) {
                    if (next_iterator[num]->next != NULL) {
                        next_iterator.emplace_back(next_iterator[num]->next);
                    }
                }
                next_iterator.erase(next_iterator.begin() + num);
                size--;
                num--;
            }
        }
        newick_order--;
    }
    gsl_matrix_free(eigenvector_temp);
    gsl_vector_free(eigenvalue_temp);
}