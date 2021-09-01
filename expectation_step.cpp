//
// Created by 박석환 on 2021/08/13.
//

#include <algorithm>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_linalg.h>

void get_eigenvector_and_inverse(gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *eigenvector, gsl_matrix *eigenvec_inverse, gsl_vector *eigenvalue) {
    //set qmatrix_temp & normalize
    gsl_matrix *qmatrix_temp = gsl_matrix_alloc(64, 64);
    gsl_matrix_memcpy(qmatrix_temp, qmatrix);
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            gsl_matrix_set(qmatrix_temp, row, col, gsl_matrix_get(qmatrix_temp, row, col) * codon_freq[col]);
        }
    }

    double sum = 0.0;
    for (size_t dia = 0; dia < 64; dia++) {
        for (size_t nondia = 0; nondia < 64; nondia++) {
            if (dia != nondia) {
                sum -= gsl_matrix_get(qmatrix_temp, dia, nondia);
            }
        }
        gsl_matrix_set(qmatrix_temp, dia, dia, sum);
        sum = 0.0;
    }

    for (size_t dia = 0; dia < 64; dia++) {
        sum -= gsl_matrix_get(qmatrix_temp, dia, dia) * codon_freq[dia];
    }
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            gsl_matrix_set(qmatrix_temp, row, col, gsl_matrix_get(qmatrix_temp, row, col)/sum);
        }
    }

    //get eigenvalues and left & right eigenvectors
    gsl_vector_complex *eigenvalue_com = gsl_vector_complex_alloc(64);//freed
    gsl_eigen_nonsymmv_workspace *eigen_space = gsl_eigen_nonsymmv_alloc(64);//freed
    gsl_matrix_complex *eigenvector_com = gsl_matrix_complex_alloc(64, 64);//freed
    gsl_eigen_nonsymmv(qmatrix_temp, eigenvalue_com, eigenvector_com, eigen_space);
    gsl_eigen_nonsymmv_free(eigen_space);
    gsl_matrix_free(qmatrix_temp);

    // check if all the eigenvalues are real && make exponential of diagonal matrix
    //Todo do I have to correct machine precision? Would it affect the result a lot? (Corrected just for case)
    gsl_complex check;
    for (size_t i = 0; i < 64; i++) {
        check = gsl_vector_complex_get(eigenvalue_com, i);
        assert(check.dat[1] == 0 && "Eigenvalues are not real");
        if (abs(check.dat[0]) < 1.0e-14) {
            gsl_vector_set(eigenvalue, i, 0);
        } else {
            gsl_vector_set(eigenvalue, i, check.dat[0]);
        }
    }

    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            gsl_matrix_set(eigenvector, row, col, gsl_matrix_complex_get(eigenvector_com, row, col).dat[0]);
        }
    }
    gsl_matrix_complex_free(eigenvector_com);

    gsl_matrix_memcpy(eigenvec_inverse, eigenvector);
    gsl_permutation *p = gsl_permutation_alloc(64);//freed
    int here = 0;
    int *signum = &here;
    gsl_linalg_LU_decomp(eigenvec_inverse, p, signum);
    gsl_linalg_LU_invx(eigenvec_inverse, p);
    gsl_permutation_free(p);
}

double cal_prob_from_qmatrix(gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, char base_before, char base_after, float branch_length) {
    for (size_t num = 0; num < 64; num++) {
        if (gsl_vector_get(eigenvalue, num) == 0) {
            continue;
        } else {
            gsl_vector_set(eigenvalue, num, pow(M_E, gsl_vector_get(eigenvalue, num) * branch_length));
        }
    }

    //compute P * e^tD * P^-1
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            gsl_matrix_set(eigenvector, row, col, gsl_matrix_get(eigenvector, row, col) * gsl_vector_get(eigenvalue, col));
        }
    }
    gsl_matrix *expon_qmatrix = gsl_matrix_alloc(64, 64);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigenvector, eigen_inverse, 0.0, expon_qmatrix);

    int row = base_before;
    int col = base_after;
    return gsl_matrix_get(expon_qmatrix, row, col);
}

double felsenstein_algorithm(newick_graph *node, char base, newick_start &start, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue) {
    if (node->previous.empty()) {
        if (node->base[base] == true) {
            node->felsenstein[base] = 1.0;
            return 1.0;
        } else {
            return 0.0;
        }
    } else {
        double sum_first = 0.0;
        double sum_second = 0.0;
        for (char num = 0; num < 64; num++) {
            sum_first += node->previous[0]->felsenstein[num] * cal_prob_from_qmatrix(eigenvector, eigen_inverse, eigenvalue, num, base, node->previous[0]->branch_length);
            sum_second += node->previous[1]->felsenstein[num] * cal_prob_from_qmatrix(eigenvector, eigen_inverse, eigenvalue, num, base, node->previous[1]->branch_length);
        }
        node->felsenstein[base] = sum_first * sum_second;
        return sum_first * sum_second;
    }
}

void conduct_felsenstein(newick_start &start) {

}