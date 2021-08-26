//
// Created by 박석환 on 2021/08/13.
//

#include <algorithm>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>


double cal_prob_from_qmatrix(gsl_matrix *qmatrix, double codon_freq, char base_before, char base_after, float branch_length) {
    //get eigenvalues and left & right eigenvectors
    //set memories
    gsl_matrix *qmatrix_temp = gsl_matrix_alloc(64, 64);
    gsl_matrix_memcpy(qmatrix_temp, qmatrix);
    gsl_vector_complex *eigenvalue = gsl_vector_complex_alloc(64);
    gsl_eigen_nonsymm_workspace *eigen_space = gsl_eigen_nonsymm_alloc(64);//freed
    gsl_matrix *Z = gsl_matrix_alloc(64, 64);

    // Schur decomposition
    gsl_eigen_nonsymm_params(1, 0, eigen_space);
    gsl_eigen_nonsymm_Z(qmatrix_temp, eigenvalue, Z, eigen_space);
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            if (row > col) {
                gsl_matrix_set(qmatrix_temp, row, col, 0.0);
            }
        }
    }

    //check if all the eigenvalues are real
    gsl_complex check;
    for (size_t i = 0; i < 64; i++) {
        check = gsl_vector_complex_get(eigenvalue, i);
        assert(check.dat[1] == 0 && "Eigenvalues are not real");
    }

    gsl_eigen_nonsymm_free(eigen_space);
    gsl_vector_complex_free(eigenvalue);

    //test diagonalization
    gsl_matrix *temp_matrix = gsl_matrix_alloc(64, 64);
    gsl_matrix_set_all(temp_matrix, 0.0);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Z, qmatrix_temp, 0.0, temp_matrix);
    gsl_matrix *result = gsl_matrix_alloc(64, 64);
    gsl_matrix_set_all(result, 0.0);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, temp_matrix, Z, 0.0, result);
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            assert(abs(gsl_matrix_get(result, row, col) - gsl_matrix_get(qmatrix, row, col)) < 1.0e-12);
        }
    }

    return 0.0;
}