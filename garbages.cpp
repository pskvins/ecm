//
// Created by 박석환 on 2021/10/13.
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
#include <limits>
#include <random>
#include <chrono>
#include <gsl/gsl_multimin.h>

struct void_to_types {
    gsl_matrix *eigenvector;
    gsl_matrix *eigenvec_inverse;
    gsl_vector *eigenvalue;
    double *codon_freq;
    newick_start *start;
    int *newick_order_max;
};

void update_upper(newick_start *start, newick_graph *end, double *codon_freq, int *newick_order_max) {
    //calculate denominator
    double denominator = 0.0;
    std::vector<newick_graph*> next_iteration = start->next;
    for (char base = 0; base < 64; base++) {
        denominator += codon_freq[base] * end->felsenstein[base];
    }

    int newick_order = *newick_order_max;
    //update upper
    while (newick_order >= 0) {
        int size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            if (next_iteration[num]->order == newick_order) {
                for (char base_parent = 0; base_parent < 64; base_parent++) {
                    for (char base_child = 0; base_child < 64; base_child++) {
                        gsl_matrix_set(next_iteration[num]->updated_upper, base_parent, base_child, next_iteration[num]->felsenstein[base_child] * next_iteration[num]->upper[base_parent] *
                                                                                                    gsl_matrix_get(next_iteration[num]->expon_matrix, base_parent, base_child) /denominator);
                    }
                }
                if (std::find(next_iteration.begin(), next_iteration.end(), next_iteration[num]->next) == next_iteration.end()) {
                    if (next_iteration[num]->next != NULL) {
                        next_iteration.emplace_back(next_iteration[num]->next);
                    }
                }
                next_iteration.erase(next_iteration.begin() + num);
                num--;
                size--;
            }
        }
        newick_order--;
    }
}

/*void calculate_expectation(newick_start *start, int *newick_order_max) {
    std::vector<newick_graph*> next_iteration = start->next;
    int newick_order = *newick_order_max;
    while (newick_order >= 0) {
        int size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            if (next_iteration[num]->order == newick_order) {
                for (char base_next = 0; base_next < 64; base_next++) {
                    double denominator = 0.0;
                    for (char base_acc = 0; base_acc < 64; base_acc++) {
                        denominator += next_iteration[num]->felsenstein[base_acc] *
                                       gsl_matrix_get(next_iteration[num]->expon_matrix, base_next, base_acc);
                    }
                    //Todo: if denominator = 0, for sure numerator is 0, so make the expectation 0 (no need to calculate) --> is this the right solution?
                    if (denominator == 0) {
                        for (char base_curr = 0; base_curr < 64; base_curr++) {
                            next_iteration[num]->expectation[64 * base_next + base_curr] = 0.0;
                        }
                        continue;
                    }
                    for (char base_curr = 0; base_curr < 64; base_curr++) {
                        next_iteration[num]->expectation[64 * base_next + base_curr] =
                                next_iteration[num]->felsenstein[base_curr] *
                                gsl_matrix_get(next_iteration[num]->expon_matrix, base_next, base_curr) *
                                next_iteration[num]->next->updated_upper[base_next] / denominator;
                    }
                }
                if (std::find(next_iteration.begin(), next_iteration.end(), next_iteration[num]->next) == next_iteration.end()) {
                    if (next_iteration[num]->next != NULL) {
                        next_iteration.emplace_back(next_iteration[num]->next);
                    }
                }
                next_iteration.erase(next_iteration.begin() + num);
                num--;
                size--;
            }
        }
        newick_order--;
    }
}*/

//copying maximization_step
void calculate_derivative(gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, newick_start *start, gsl_vector *gradient, int *newick_order_max) {
    std::cout << "start calculating gradient" << std::endl;
    gsl_vector_set_all(gradient, 0.0);
    gsl_matrix *dxepon_matrix[64 * 63 / 2];
    for (int num = 0; num < 2016; num++) {
        dxepon_matrix[num] = gsl_matrix_alloc(64, 64);//freed
    }
    double diagonal[64] = {0};
    double normalize = 0.0;
    for (int row = 0; row < 64; row++) {
        for (int col = 0; col < 64; col++) {
            if (row != col) {
                diagonal[row] -= gsl_matrix_get(qmatrix, row, col) * codon_freq[col];
            }
        }
        normalize -= diagonal[row] * codon_freq[row];
    }
    thread_local int newick_order = *newick_order_max;
    thread_local std::vector<newick_graph*> next_iteration = start->next;
    thread_local size_t row_tar = 0;
    thread_local size_t col_tar = 0;
    omp_set_num_threads(4);
    //#pragma omp parallel for schedule(dynamic)
    for (size_t diff = 0; diff < 64 * 63 / 2; diff++) {
        gsl_matrix *F = gsl_matrix_alloc(64, 64);//freed
        newick_order = *newick_order_max;
        next_iteration = start->next;
        row_tar = num_to_coordinate[diff][0];
        col_tar = num_to_coordinate[diff][1];
        while (newick_order >= 0) {
            size_t size = next_iteration.size();
            for (size_t num = 0; num < size; num++) {
                if (next_iteration[num]->order == newick_order) {
                    gsl_matrix_set_all(dxepon_matrix[diff], 0.0);
                    for (size_t row = 0; row < 64; row++) {
                        for (size_t col = 0; col < 64; col++) {
                            if (row == col && (row == row_tar || row == col_tar)) {
                                gsl_matrix_set(dxepon_matrix[diff], row, col,
                                               (-2 * codon_freq[row_tar] * codon_freq[col_tar] * diagonal[row] -
                                                codon_freq[col_tar] * normalize) / pow(normalize, 2));
                            } else if ((row == row_tar && col == col_tar) || (row == col_tar && col == row_tar)) {
                                gsl_matrix_set(dxepon_matrix[diff], row, col,
                                               (codon_freq[col] * normalize -
                                                2 * codon_freq[row] * pow(codon_freq[col], 2) *
                                                gsl_matrix_get(qmatrix, row_tar, col_tar) / pow(normalize, 2)));
                            } else {
                                gsl_matrix_set(dxepon_matrix[diff], row, col,
                                               -gsl_matrix_get(qmatrix, row, col) * 2 * codon_freq[row_tar] *
                                               codon_freq[col_tar] / pow(normalize, 2));
                            }
                        }
                    }
                    gsl_matrix *temp_matrix = gsl_matrix_alloc(64, 64);//freed
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigen_inverse, dxepon_matrix[diff], 0.0, temp_matrix);
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, next_iteration[num]->branch_length, temp_matrix, eigenvector, 0.0, dxepon_matrix[diff]);
                    for (size_t row = 0; row < 64; row++) {
                        for (size_t col = 0; col < 64; col++) {
                            if (abs(gsl_vector_get(eigenvalue, row) - gsl_vector_get(eigenvalue, col)) < std::numeric_limits<double>::epsilon()) {
                                gsl_matrix_set(F, row, col,
                                               next_iteration[num]->branch_length * pow(M_E, gsl_vector_get(eigenvalue, row) *
                                                                                             next_iteration[num]->branch_length));
                            } else {
                                gsl_matrix_set(F, row, col,
                                               (pow(M_E, gsl_vector_get(eigenvalue, row) *
                                                         next_iteration[num]->branch_length)
                                                - pow(M_E, gsl_vector_get(eigenvalue, col) *
                                                           next_iteration[num]->branch_length))
                                               / (gsl_vector_get(eigenvalue, row) - gsl_vector_get(eigenvalue, col)));
                            }
                        }
                    }
                    gsl_matrix_mul_elements(dxepon_matrix[diff], F);
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigenvector, dxepon_matrix[diff], 0.0, temp_matrix);
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp_matrix, eigen_inverse, 0.0, dxepon_matrix[num]);
                    gsl_matrix_free(temp_matrix);
                    //putting negative of gradient
                    for (size_t row = 0; row < 64; row++) {
                        for (size_t col = 0; col < 64; col++) {
                            gsl_vector_set(gradient, diff, gsl_vector_get(gradient, diff) - next_iteration[num]->expectation[64 * row + col]
                                                                                            * gsl_matrix_get(dxepon_matrix[diff], row, col) / gsl_matrix_get(next_iteration[num]->expon_matrix[0], row, col));
                            if (gsl_matrix_get(next_iteration[num]->expon_matrix[0], row, col) == 0) {
                                std::cout << "cannot divide by zero!" << std::endl;
                                exit(-1);
                            }
                        }
                    }
                    if (std::find(next_iteration.begin(), next_iteration.end(), next_iteration[num]->next) == next_iteration.end()) {
                        if (next_iteration[num]->next != NULL) {
                            next_iteration.emplace_back(next_iteration[num]->next);
                        }
                    }
                    next_iteration.erase(next_iteration.begin() + num);
                    num--;
                    size--;
                }
            }
            newick_order--;
        }
        gsl_matrix_free(F);
    }
    for (int num = 0; num < 2016; num++) {
        gsl_matrix_free(dxepon_matrix[num]);
    }
    for (int i = 0; i < 2016; i++) {
        if (isnan(gsl_vector_get(gradient, i))) {
            std::cout << "gradient is nan" << std::endl;
            exit(-1);
        }
    }
    std::cout << "end calculating gradient" << std::endl;
}

void update_and_normalize_qmatirx(gsl_matrix *qmatrix, double *codon_freq, gsl_vector *direction, double lambda, gsl_matrix *new_matrix, bool &non_diag_neg) {
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            if (row > col) {
                gsl_matrix_set(new_matrix, row, col, (gsl_matrix_get(qmatrix, row, col) + lambda * gsl_vector_get(direction,(row - 1) * (row) / 2 + col)) * codon_freq[col]);
                gsl_matrix_set(new_matrix, col, row, (gsl_matrix_get(qmatrix, col, row) + lambda * gsl_vector_get(direction,(col - 1) * (col) / 2 + row)) * codon_freq[row]);
            } else if (row == col) {
                gsl_matrix_set(new_matrix, row, col, 0.0);
            }
        }
    }
    //check if new matrix's non-diagonal is negative
    if (gsl_matrix_isnonneg(new_matrix) == 0) {
        non_diag_neg = true;
        std::cout << "have to start again" << std::endl;
    }
    //if there is negative entry for non-diagonal, reset the qmatrix to be the fixed version of new_matrix
    if (non_diag_neg == true) {
        for (size_t row = 0; row < 64; row++) {
            for (size_t col = 0; col < 64; col++) {
                if (row > col) {
                    gsl_matrix_set(qmatrix, row, col, gsl_matrix_get(qmatrix, row, col) + lambda * gsl_vector_get(direction, (row - 1) * (row) / 2 + col));
                    gsl_matrix_set(qmatrix, col, row, gsl_matrix_get(qmatrix, row, col));
                }
            }
        }
        double sum = 0.0;
        for (size_t dia = 0; dia < 64; dia++) {
            for (size_t nondia = 0; nondia < 64; nondia++) {
                if (dia != nondia) {
                    sum -= gsl_matrix_get(qmatrix, dia, nondia);
                }
            }
            gsl_matrix_set(qmatrix, dia, dia, sum);
            sum = 0.0;
        }
        changing_to_nonneg_matrix(qmatrix);
        return;
    }
    double sum = 0.0;
    for (size_t dia = 0; dia < 64; dia++) {
        for (size_t nondia = 0; nondia < 64; nondia++) {
            if (dia != nondia) {
                sum -= gsl_matrix_get(new_matrix, dia, nondia);
            }
        }
        gsl_matrix_set(new_matrix, dia, dia, sum);
        sum = 0.0;
    }

    for (size_t dia = 0; dia < 64; dia++) {
        sum -= gsl_matrix_get(new_matrix, dia, dia) * codon_freq[dia];
    }
    gsl_matrix_scale(new_matrix, 1 / sum);
}

double calculate_maximization_function(newick_start *start, short int index, int *newick_order_max) {
    double result = 0;
    std::vector<newick_graph*> next_iteration = start->next;
    int newick_order = *newick_order_max;
    while(newick_order >= 0) {
        size_t size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            if (next_iteration[num]->order == newick_order) {
                for (size_t row = 0; row < 64; row++) {
                    for (size_t col = 0; col < 64; col++) {
                        result -= next_iteration[num]->expectation[64 * row + col] *
                                  log(gsl_matrix_get(next_iteration[num]->expon_matrix[index], row, col));
                        if (gsl_matrix_get(next_iteration[num]->expon_matrix[index], row, col) <= 0) {
                            std::cout << "log0 or log(neg)!! : " << gsl_matrix_get(next_iteration[num]->expon_matrix[index], row, col) << "," << log(gsl_matrix_get(next_iteration[num]->expon_matrix[index], row, col)) << std::endl;
                            //exit(-1);
                        }
                    }
                }
                if (std::find(next_iteration.begin(), next_iteration.end(), next_iteration[num]->next) == next_iteration.end()) {
                    if (next_iteration[num]->next != NULL) {
                        next_iteration.emplace_back(next_iteration[num]->next);
                    }
                }
                next_iteration.erase(next_iteration.begin() + num);
                num--;
                size--;
            }
        }
        newick_order--;
    }
    if (isnan(result)) {
        std::cout << "function value is nan" << std::endl;
        exit(-1);
    }
    std::cout << "function value : " << result << std::endl;
    return result;
}

void update_expon_matrix(newick_start *start, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, int *newick_order_max) {
    std::vector<newick_graph*> next_iterator = start->next;
    gsl_matrix *eigenvector_temp = gsl_matrix_alloc(64, 64);//freed
    gsl_vector *eigenvalue_temp = gsl_vector_alloc(64);//freed
    int newick_order = *newick_order_max;
    while (newick_order >= 0) {
        size_t size = next_iterator.size();
        for (size_t num = 0; num < size; num++) {
            if (next_iterator[num]->order == newick_order) {
                gsl_matrix_memcpy(next_iterator[num]->expon_matrix[1],
                                  calculate_expon_matrix(eigenvector, eigen_inverse, eigenvalue, eigenvector_temp, eigenvalue_temp, next_iterator[num]->branch_length));
                if (std::find(next_iterator.begin(), next_iterator.end(), next_iterator[num]->next) == next_iterator.end()) {
                    if (next_iterator[num]->next != NULL) {
                        next_iterator.emplace_back(next_iterator[num]->next);
                    }
                }
                next_iterator.erase(next_iterator.begin() + num);
                num--;
                size--;
            }
        }
        newick_order--;
    }
    gsl_matrix_free(eigenvector_temp);
    gsl_vector_free(eigenvalue_temp);
}

/*void eigendecomp_for_backtracking(gsl_matrix *qmatrix_temp, gsl_matrix *eigenvector, gsl_matrix *eigenvec_inverse, gsl_vector *eigenvalue,
                                  double *codon_freq, gsl_matrix *qmatrix, double &lamb, double lambmin, gsl_vector *direction, bool &lamb_smaller_min) {
    //get eigenvalues and left & right eigenvectors
    gsl_vector_complex *eigenvalue_com = gsl_vector_complex_alloc(64);//freed
    gsl_eigen_nonsymmv_workspace *eigen_space = gsl_eigen_nonsymmv_alloc(64);//freed
    gsl_matrix_complex *eigenvector_com = gsl_matrix_complex_alloc(64, 64);//freed
    bool status_non_zero = false;
    gsl_set_error_handler_off();
    int status = gsl_eigen_nonsymmv(qmatrix_temp, eigenvalue_com, eigenvector_com, eigen_space);
    std::cout << status << std::endl;
    if (status == 11) {
        status_non_zero = true;
        while (status != 0 && lamb > lambmin) {
            lamb *= 0.9;
            update_and_normalize_qmatirx(qmatrix, codon_freq, direction, lamb, qmatrix_temp);
            status = gsl_eigen_nonsymmv(qmatrix_temp, eigenvalue_com, eigenvector_com, eigen_space);
        }
        if (lamb <= lambmin && status_non_zero == true) {
            lamb_smaller_min = true;
            return;
        }
    } else if (status != 0) {
        std::cout << "status not 0! : " << status << std::endl;
    }
    gsl_set_error_handler(NULL);

    // check if all the eigenvalues are real && make exponential of diagonal matrix
    //Todo: do I have to correct machine precision? Would it affect the result a lot? (Corrected just in case)
    gsl_complex check;
    bool check_real = true;
    for (size_t i = 0; i < 64; i++) {
        check = gsl_vector_complex_get(eigenvalue_com, i);
        if (check.dat[1] != 0 ) {
            std::cout << "Eigenvalues are not real : " << check.dat[1] << std::endl;
            check_real = false;
            update_and_normalize_qmatirx(qmatrix, codon_freq, direction, lamb, qmatrix_temp);
            for (size_t row = 0; row < 64; row++) {
                for (size_t col = 0; col < 64; col++) {
                    if (row > col) {
                        std::cout << gsl_matrix_get(qmatrix_temp, row, col) << ", ";
                    }
                }
                std::cout << std::endl;
            }
            throw("END!");
            break;
        }
        if (abs(check.dat[0]) < 1.0e-14) {
            gsl_vector_set(eigenvalue, i, 0);
        } else {
            gsl_vector_set(eigenvalue, i, check.dat[0]);
        }
    }
    while (check_real == false && lamb > lambmin) {
        lamb *= 0.9;
        update_and_normalize_qmatirx(qmatrix, codon_freq, direction, lamb, qmatrix_temp);
        status = gsl_eigen_nonsymmv(qmatrix_temp, eigenvalue_com, eigenvector_com, eigen_space);
        if (status == 11) {
            while (status != 0 && lamb > lambmin) {
                lamb *= 0.9;
                update_and_normalize_qmatirx(qmatrix, codon_freq, direction, lamb, qmatrix_temp);
                status = gsl_eigen_nonsymmv(qmatrix_temp, eigenvalue_com, eigenvector_com, eigen_space);
            }
        } else if (status != 0) {
            std::cout << "status not 0! : " << status << std::endl;
        }
        for (size_t i = 0; i < 64; i++) {
            check = gsl_vector_complex_get(eigenvalue_com, i);
            if (check.dat[1] != 0 ) {
                std::cout << "Eigenvalues are not real" << std::endl;
                check_real = false;
                break;
            }
            if (abs(check.dat[0]) < 1.0e-14) {
                gsl_vector_set(eigenvalue, i, 0);
            } else {
                gsl_vector_set(eigenvalue, i, check.dat[0]);
            }
            if (i == 63) {
                check_real = true;
            }
        }
    }
    if (lamb <= lambmin && status_non_zero == true) {
        lamb_smaller_min = true;
        return;
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
    gsl_eigen_nonsymmv_free(eigen_space);
    gsl_vector_complex_free(eigenvalue_com);
}*/


//Todo: do I have to limit the length of each step? If I do, then how much? -> first try : 60 (too big result), second try : 100 (got bigger)
double backtracking(newick_start *start, gsl_matrix *qmatrix, double *codon_freq, gsl_vector *direction, gsl_vector *gradient, double &lamb,
                    gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, int *newick_order_max, bool &non_diag_neg) {
    gsl_matrix *qmatrix_temp = gsl_matrix_alloc(64, 64);//freed
    gsl_vector *gradient_new = gsl_vector_alloc(2016);//freed
    //TODO: Do I have to normalize qmatrix_new in this step??-->seems like has to be normalized (equation needs normalization)
    double sum = 0.0;
    for (int i = 0; i < 2016; i++) {
        sum += pow(gsl_vector_get(direction, i), 2);
    }
    double direction_max = 0;
    for (int i = 0; i < 2016; i++) {
        if (direction_max < abs(gsl_vector_get(direction, i))) {
            direction_max = abs(gsl_vector_get(direction, i));
        }
    }
    /*sum = sqrt(sum);
    if (sum > 2016) {
        gsl_vector_scale(direction, 2016 / sum);
    }*/
    gsl_vector_scale(direction, 1 / (10 * direction_max));
    double funcg_ori = calculate_maximization_function(start, 0, newick_order_max);
    std::cout << "original function value : " << funcg_ori << std::endl;
    double slope = 0.0;
    for (int num = 0; num < 2016; num++) {
        slope += gsl_vector_get(gradient, num) * gsl_vector_get(direction, num);
    }
    if (slope >= 0.0) {
        std::cout << "Direction not going down, slope : " << slope << std::endl;
    } else {
        std::cout << "Direction going down, slope : " << slope << std::endl;
    }
    double pseudo_slope = 0.0;
    const double machine_precision = std::numeric_limits<double>::epsilon();
    const double alpha = 1.0e-4;
    double temp;
    double test = 0.0;
    size_t row;
    size_t col;
    for (size_t num = 0; num < 63 * 64 / 2; num++) {
        row = num_to_coordinate[num][0];
        col = num_to_coordinate[num][1];
        temp = abs(gsl_vector_get(direction, num)) / std::max(abs(gsl_matrix_get(qmatrix, row, col)), 1.0);
        if (temp > test) {
            test = temp;
        }
    }
    double lambmin = machine_precision / test;
    double lamb_old = 0.0;
    double lamb_tmp, cubic_a, cubic_b, cubic_comp_first, cubic_comp_secnd;
    double funcg;
    double funcg_old = 0.0;
    double inside_sqrt;
    while (true) {
        std::cout << "current lambda : " << lamb << std::endl;
        //normalize after adding direction
        update_and_normalize_qmatirx(qmatrix, codon_freq, direction, lamb, qmatrix_temp, non_diag_neg);
        if (non_diag_neg == true) {
            return NAN;
        }

        bool lamb_smaller_min = false;
        //get eigenvalues and eigenvectors
        get_eigenvector_and_inverse(qmatrix_temp, eigenvector, eigen_inverse, eigenvalue);
        if (lamb_smaller_min == true) {
            std::cout << "so sad situation, lamb : " << lamb << std::endl;
            exit (-1);
        }

        //update matrix in tree
        update_expon_matrix(start, eigenvector, eigen_inverse, eigenvalue, newick_order_max);

        //calculate g(lambda)
        funcg = calculate_maximization_function(start, 1, newick_order_max);

        if (lamb < lambmin) {
            gsl_matrix_free(qmatrix_temp);
            lamb = 0.0;
            std::cout << "lamb smaller than min : " << lamb << std::endl;
            gsl_vector_free(gradient_new);
            return funcg;
        } else if (funcg <= funcg_ori + alpha * lamb * slope) {
            //calculate new derivative
            for (size_t row = 0; row < 64; row++) {
                for (size_t col = 0; col < 64; col++) {
                    if (row > col) {
                        gsl_matrix_set(qmatrix_temp, row, col, (gsl_matrix_get(qmatrix, row, col) + lamb * gsl_vector_get(direction,(row - 1) * (row) / 2 + col)) * codon_freq[col]);
                        gsl_matrix_set(qmatrix_temp, col, row, (gsl_matrix_get(qmatrix, col, row) + lamb * gsl_vector_get(direction,(col - 1) * (col) / 2 + row)) * codon_freq[row]);
                    }
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
            calculate_derivative(qmatrix_temp, codon_freq, eigenvector, eigen_inverse, eigenvalue, start, gradient_new, newick_order_max);
            //get pseudo_slope
            for (int num = 0; num < 2016; num++) {
                pseudo_slope += gsl_vector_get(gradient_new, num) * gsl_vector_get(direction, num);
            }
            if (abs(pseudo_slope) > 0.9 * abs(slope)) {
                if (lamb == 1.0) {
                    lamb_tmp = -slope / (2.0 * (funcg - funcg_ori - slope));
                } else {
                    cubic_comp_first = funcg - funcg_ori - lamb * slope;
                    cubic_comp_secnd = funcg_old - funcg_ori - lamb_old * slope;
                    cubic_a = (cubic_comp_first / (lamb * lamb) - cubic_comp_secnd / (lamb_old * lamb_old)) / (lamb - lamb_old);
                    cubic_b = (-lamb_old * cubic_comp_first / (lamb * lamb) + lamb * cubic_comp_secnd / (lamb_old * lamb_old) / (lamb - lamb_old));
                    if (cubic_a == 0.0) {
                        lamb_tmp = -slope / (2.0 * cubic_b);
                    } else {
                        inside_sqrt = cubic_b * cubic_b - 3.0 * cubic_a * slope;
                        //Todo: what if inside_sqrt goes to positive infinite? Is it a possible situation?
                        if (inside_sqrt < 0.0 || isnan(inside_sqrt)) {
                            lamb_tmp = 0.5 * lamb;
                        } else if (cubic_b <= 0.0) {
                            lamb_tmp = (-cubic_b + sqrt(inside_sqrt)) / (3.0 * cubic_a);
                        } else {
                            lamb_tmp = -slope / (cubic_b + sqrt(inside_sqrt));
                        }
                    }
                    if (lamb_tmp > 0.5 * lamb) {
                        lamb_tmp = 0.5 * lamb;
                    }
                }
                lamb_old = lamb;
                funcg_old = funcg;
                lamb = std::max(lamb_tmp, 0.1 * lamb);
                continue;
            }
            gsl_matrix_free(qmatrix_temp);
            gsl_vector_free(gradient_new);
            std::cout << "lamb : " << lamb << std::endl;
            return funcg;
        } else {
            if (lamb == 1.0) {
                lamb_tmp = -slope / (2.0 * (funcg - funcg_ori - slope));
            } else {
                cubic_comp_first = funcg - funcg_ori - lamb * slope;
                cubic_comp_secnd = funcg_old - funcg_ori - lamb_old * slope;
                cubic_a = (cubic_comp_first / (lamb * lamb) - cubic_comp_secnd / (lamb_old * lamb_old)) / (lamb - lamb_old);
                cubic_b = (-lamb_old * cubic_comp_first / (lamb * lamb) + lamb * cubic_comp_secnd / (lamb_old * lamb_old) / (lamb - lamb_old));
                if (cubic_a == 0.0) {
                    lamb_tmp = -slope / (2.0 * cubic_b);
                } else {
                    inside_sqrt = cubic_b * cubic_b - 3.0 * cubic_a * slope;
                    //Todo: what if inside_sqrt goes to positive infinite? Is it a possible situation?
                    if (inside_sqrt < 0.0 || isnan(inside_sqrt)) {
                        lamb_tmp = 0.5 * lamb;
                    } else if (cubic_b <= 0.0) {
                        lamb_tmp = (-cubic_b + sqrt(inside_sqrt)) / (3.0 * cubic_a);
                    } else {
                        lamb_tmp = -slope / (cubic_b + sqrt(inside_sqrt));
                    }
                }
                if (lamb_tmp > 0.5 * lamb) {
                    lamb_tmp = 0.5 * lamb;
                }
            }
        }
        lamb_old = lamb;
        funcg_old = funcg;
        lamb = std::max(lamb_tmp, 0.1 * lamb);
    }
}

double quasi_Newton_method(newick_start *start, gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, int *newick_order_max, bool &non_diag_neg) {
    std::cout << "starting maximization step" << std::endl;
    const double machine_precision = 2.0e-14;
    const int iterate_max = 200;
    //initialize matrix
    gsl_matrix *hessian = gsl_matrix_alloc(2016, 2016);//freed
    gsl_matrix *hessian_append = gsl_matrix_alloc(2016, 2016);//freed
    gsl_matrix *hessian_temp = gsl_matrix_alloc(2016, 2016);//freed
    gsl_matrix *hessian_vector = gsl_matrix_alloc(1, 2016);//freed
    gsl_matrix *yk = gsl_matrix_alloc(2016, 1);//freed
    gsl_matrix_set_all(hessian, 0.0);
    for (int num = 0; num < 2016; num++) {
        gsl_matrix_set(hessian, num, num, 1.0);
    }
    gsl_vector *gradient = gsl_vector_alloc(2016);//freed
    gsl_vector *gradient_new = gsl_vector_alloc(2016);//freed
    gsl_vector_set_all(gradient, 0.0);
    gsl_vector_set_all(gradient_new, 0.0);
    //get gradient
    calculate_derivative(qmatrix, codon_freq, eigenvector, eigen_inverse, eigenvalue, start, gradient, newick_order_max);

    gsl_vector *dx = gsl_vector_alloc(2016);//freed
    gsl_matrix *x = gsl_matrix_alloc(2016, 1);//freed
    gsl_matrix *x_new = gsl_matrix_alloc(2016, 1);//freed
    //delta x and x for first iteration
    for (int row = 0; row < 2016; row++) {
        gsl_vector_set(dx, row, -gsl_vector_get(gradient, row));
    }
    for (int row = 0; row < 64; row++) {
        for (int col = 0; col < 64; col++) {
            if (row > col) {
                gsl_matrix_set(x, (row - 1) * (row) / 2 + col, 0, gsl_matrix_get(qmatrix, row, col));
            }
        }
    }

    int row_map, col_map;
    gsl_matrix *qmatrix_new = gsl_matrix_alloc(64, 64);//freed
    double function_value;
    //start loop
    for (int its = 0; its < iterate_max; its++) {
        std::cout << std::endl;
        std::cout << its << "th iteration" << std::endl;
        double lambda = 1.0;
        function_value = backtracking(start, qmatrix, codon_freq, dx, gradient, lambda, eigenvector, eigen_inverse, eigenvalue, newick_order_max, non_diag_neg);
        if (non_diag_neg == true) {
            std::cout << "end maximization step" << std::endl;
            return NAN;
        }
        for (int num = 0; num < 2016; num++) {
            gsl_vector_set(dx, num, lambda * gsl_vector_get(dx, num));
            row_map = num_to_coordinate[num][0];
            col_map = num_to_coordinate[num][1];
            gsl_matrix_set(qmatrix_new, row_map, col_map, gsl_matrix_get(qmatrix, row_map, col_map) + gsl_vector_get(dx, num));
            gsl_matrix_set(qmatrix_new, col_map, row_map, gsl_matrix_get(qmatrix_new, row_map, col_map));
            gsl_matrix_set(x_new, num, 0, gsl_matrix_get(x, num, 0) + gsl_vector_get(dx, num));
        }
        //set rest of qmatrix_new and x_new(already set)
        double sum = 0.0;
        for (size_t dia = 0; dia < 64; dia++) {
            for (size_t nondia = 0; nondia < 64; nondia++) {
                if (dia != nondia) {
                    sum -= gsl_matrix_get(qmatrix_new, dia, nondia);
                }
            }
            gsl_matrix_set(qmatrix_new, dia, dia, sum);
            sum = 0.0;
        }

        double test = 0.0;
        for (int num = 0; num < 2016; num++) {
            double temp = abs(gsl_vector_get(dx, num)) / std::max(abs(gsl_matrix_get(x_new, num, 0)), 1.0);
            if (temp > test) {
                test = temp;
            }
        }
        if (test < 4 * machine_precision) {
            gsl_matrix_memcpy(qmatrix, qmatrix_new);
            gsl_matrix_free(hessian);
            gsl_matrix_free(hessian_append);
            gsl_matrix_free(hessian_temp);
            gsl_matrix_free(hessian_vector);
            gsl_vector_free(gradient);
            gsl_vector_free(gradient_new);
            gsl_matrix_free(x);
            gsl_matrix_free(x_new);
            gsl_vector_free(dx);
            gsl_matrix_free(qmatrix_new);
            gsl_matrix_free(yk);
            //change all the negative entries into non-negative entries (Israel, Rosenthal & Wei (2001))
            changing_to_nonneg_matrix(qmatrix);
            std::cout << "end maximization step" << std::endl;
            return function_value;
        }
        calculate_derivative(qmatrix_new, codon_freq, eigenvector, eigen_inverse, eigenvalue, start, gradient_new, newick_order_max);
        test = 0.0;
        for (int num = 0; num < 2016; num++) {
            double temp = abs(gsl_vector_get(gradient_new, num)) / std::max(abs(gsl_matrix_get(x_new, num, 0)), 1.0);
            if (temp > test) {
                test = temp;
            }
        }
        if (test < machine_precision) {
            gsl_matrix_memcpy(qmatrix, qmatrix_new);
            gsl_matrix_free(hessian);
            gsl_matrix_free(hessian_append);
            gsl_matrix_free(hessian_temp);
            gsl_matrix_free(hessian_vector);
            gsl_vector_free(gradient);
            gsl_vector_free(gradient_new);
            gsl_matrix_free(x);
            gsl_matrix_free(x_new);
            gsl_vector_free(dx);
            gsl_matrix_free(qmatrix_new);
            gsl_matrix_free(yk);
            //change all the negative entries into non-negative entries (Israel, Rosenthal & Wei (2001))
            changing_to_nonneg_matrix(qmatrix);
            std::cout << "end maximization step" << std::endl;
            return function_value;
        }
        //now calculate yk
        for (int num = 0 ; num < 2016; num++) {
            gsl_matrix_set(yk, num, 0, gsl_vector_get(gradient_new, num) - gsl_vector_get(gradient, num));
        }
        //update hessian matrix
        double denominator = 0.0;
        for (int num = 0; num < 2016; num++) {
            denominator += gsl_matrix_get(yk, num, 0) * gsl_vector_get(dx, num);
        }
        std::cout << "Denominator value : " << denominator << std::endl;
        double scalar_bfgs = 0.0;
        double scalar_temp = 0.0;
        for (int num = 0; num < 2016 ; num++) {
            for (int right = 0; right < 2016; right++) {
                scalar_temp += gsl_matrix_get(hessian, num, right) * gsl_matrix_get(yk, right, 0);
            }
            scalar_bfgs += scalar_temp * gsl_matrix_get(yk, num, 0);
            scalar_temp = 0.0;
        }
        scalar_bfgs += denominator;
        gsl_matrix_set_all(hessian_append, 0.0);
        gsl_blas_dsyr(CblasLower, scalar_bfgs / (denominator * denominator), dx, hessian_append);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, yk, hessian, 0.0, hessian_vector);
        gsl_matrix_add(hessian, hessian_append);
        gsl_matrix *dx_temp = gsl_matrix_alloc(2016, 1);//freed
        for (int i = 0; i < 2016; i++) {
            gsl_matrix_set(dx_temp, i, 0, gsl_vector_get(dx, i));
        }
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1 / denominator, dx_temp, hessian_vector, 0.0, hessian_append);
        gsl_matrix_free(dx_temp);
        gsl_matrix_transpose_memcpy(hessian_temp, hessian_append);
        gsl_matrix_add(hessian_append, hessian_temp);
        gsl_matrix_add(hessian, hessian_append);

        //old version of updating, unable it when using upper one
        /*gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0 / denominator, dx, gradient, 0.0, hessian_append);
        for (int num = 0; num < 2016; num++) {
            gsl_matrix_set(hessian_append, num, num, 1.0 + gsl_matrix_get(hessian_append, num, num));
        }
        gsl_matrix *hessian_temp = gsl_matrix_alloc(2016, 2016);
        gsl_matrix_memcpy(hessian_temp, hessian);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, hessian_append, hessian_temp, 0.0, hessian);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0 / denominator, gradient, dx, 0.0, hessian_append);
        for (int num = 0; num < 2016; num++) {
            gsl_matrix_set(hessian_append, num, num, 1.0 + gsl_matrix_get(hessian_append, num, num));
        }
        gsl_matrix_memcpy(hessian_temp, hessian);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, hessian_temp, hessian_append, 0.0, hessian);
        gsl_matrix_free(hessian_temp);
        for (int row = 0; row < 2016; row++) {
            for (int col = 0; col < 2016; col++) {
                if (row > col) {
                    gsl_matrix_set(hessian_append, row, col,gsl_matrix_get(dx, row, 0) * gsl_matrix_get(dx, col, 0) / denominator);
                    gsl_matrix_set(hessian_append, col, row, gsl_matrix_get(hessian_append, row, col));
                } else if (row == col) {
                    gsl_matrix_set(hessian_append, row, col,gsl_matrix_get(dx, row, 0) * gsl_matrix_get(dx, col, 0) / denominator);
                }
            }
        }
        gsl_matrix_add(hessian, hessian_append);*/

        //update gradient, qmatirx x and dx
        gsl_matrix_memcpy(x, x_new);
        gsl_vector_memcpy(gradient, gradient_new);
        gsl_matrix_memcpy(qmatrix, qmatrix_new);
        gsl_blas_dgemv(CblasNoTrans, -1.0, hessian, gradient, 0.0, dx);
    }
    std::cout << "Too many iterations" << std::endl;
    exit(-1);
}

//functions needed for gsl quasi newton method
double get_function_value_max_step(const gsl_vector *x, void* params) {
    //allocate x and normalize the matrix
    void_to_types* params_real = (void_to_types*) params;
    gsl_matrix *qmatrix = gsl_matrix_alloc(64, 64);//freed
    for (int row = 0; row < 64; row++) {
        for (int col = 0; col < 64; col++) {
            if (row > col) {
                gsl_matrix_set(qmatrix, row, col, gsl_vector_get(x,(row - 1) * (row) / 2 + col) * params_real->codon_freq[col]);
                gsl_matrix_set(qmatrix, col, row, gsl_vector_get(x,(row - 1) * (row) / 2 + col) * params_real->codon_freq[row]);
            }
        }
    }
    double sum = 0.0;
    for (size_t dia = 0; dia < 64; dia++) {
        for (size_t nondia = 0; nondia < 64; nondia++) {
            if (dia != nondia) {
                sum -= gsl_matrix_get(qmatrix, dia, nondia);
            }
        }
        gsl_matrix_set(qmatrix, dia, dia, sum);
        sum = 0.0;
    }

    for (size_t dia = 0; dia < 64; dia++) {
        sum -= gsl_matrix_get(qmatrix, dia, dia) * params_real->codon_freq[dia];
    }
    gsl_matrix_scale(qmatrix, 1 / sum);

    get_eigenvector_and_inverse(qmatrix, params_real->eigenvector, params_real->eigenvec_inverse, params_real->eigenvalue);
    gsl_matrix_free(qmatrix);

    update_expon_matrix(params_real->start, params_real->eigenvector, params_real->eigenvec_inverse, params_real->eigenvalue, params_real->newick_order_max);

    return calculate_maximization_function(params_real->start, 1, params_real->newick_order_max);
}

void get_gradient_max_step(const gsl_vector *x, void* params, gsl_vector *gradient) {
    void_to_types* params_real = (void_to_types*) params;
    gsl_matrix *qmatrix = gsl_matrix_alloc(64, 64);//freed
    for (int row = 0; row < 64; row++) {
        for (int col = 0; col < 64; col++) {
            if (row > col) {
                gsl_matrix_set(qmatrix, row, col, gsl_vector_get(x,(row - 1) * (row) / 2 + col) * params_real->codon_freq[col]);
                gsl_matrix_set(qmatrix, col, row, gsl_vector_get(x,(row - 1) * (row) / 2 + col) * params_real->codon_freq[row]);
            }
        }
    }
    double sum = 0.0;
    for (size_t dia = 0; dia < 64; dia++) {
        for (size_t nondia = 0; nondia < 64; nondia++) {
            if (dia != nondia) {
                sum -= gsl_matrix_get(qmatrix, dia, nondia);
            }
        }
        gsl_matrix_set(qmatrix, dia, dia, sum);
        sum = 0.0;
    }

    for (size_t dia = 0; dia < 64; dia++) {
        sum -= gsl_matrix_get(qmatrix, dia, dia) * params_real->codon_freq[dia];
    }
    gsl_matrix_scale(qmatrix, 1 / sum);
    get_eigenvector_and_inverse(qmatrix, params_real->eigenvector, params_real->eigenvec_inverse, params_real->eigenvalue);

    for (int row = 0; row < 64; row++) {
        for (int col = 0; col < 64; col++) {
            if (row > col) {
                gsl_matrix_set(qmatrix, row, col, gsl_vector_get(x,(row - 1) * (row) / 2 + col));
                gsl_matrix_set(qmatrix, col, row, gsl_matrix_get(qmatrix, row, col));
            }
        }
    }
    calculate_derivative(qmatrix, params_real->codon_freq, params_real->eigenvector, params_real->eigenvec_inverse, params_real->eigenvalue, params_real->start, gradient, params_real->newick_order_max);
    gsl_matrix_free(qmatrix);
}

void get_func_value_and_gradient_max_step(const gsl_vector *x, void* params, double *f, gsl_vector *gradient) {
    *f = get_function_value_max_step(x, params);
    get_gradient_max_step(x, params, gradient);
}

double quasi_newton_method_gsl(gsl_vector *x, void_to_types *params) {
    std::cout << "starting maximization step" << std::endl;
    gsl_multimin_fdfminimizer *quasi_newton_bfgs = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, 2016);//freed
    gsl_multimin_function_fdf functions;
    functions.n = 63 * 64 / 2;
    functions.f = &get_function_value_max_step;
    functions.df = &get_gradient_max_step;
    functions.fdf = &get_func_value_and_gradient_max_step;
    functions.params = (void*)params;
    gsl_multimin_fdfminimizer_set(quasi_newton_bfgs, &functions, x, 3.0, 0.9);
    int iter = 0;
    int status;
    bool nondiag_neg = false;
    do {
        iter++;
        std::cout << iter << std::endl;
        status = gsl_multimin_fdfminimizer_iterate(quasi_newton_bfgs);
        if (status == GSL_ENOPROG) {
            std::cout << "Not making progress" << std::endl;
            std::cout << "x" << std::endl;
            for (int row = 0; row < 64; row++) {
                for (int col = 0; col < 64; col++) {
                    if (row > col) {
                        std::cout << gsl_vector_get(gsl_multimin_fdfminimizer_x(quasi_newton_bfgs), (row - 1) * (row) / 2 + col) << ", ";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << "dx" << std::endl;
            for (int row = 0; row < 64; row++) {
                for (int col = 0; col < 64; col++) {
                    if (row > col) {
                        std::cout << gsl_vector_get(gsl_multimin_fdfminimizer_dx(quasi_newton_bfgs), (row - 1) * (row) / 2 + col) << ", ";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << "gradient" << std::endl;
            for (int row = 0; row < 64; row++) {
                for (int col = 0; col < 64; col++) {
                    if (row > col) {
                        std::cout << gsl_vector_get(gsl_multimin_fdfminimizer_gradient(quasi_newton_bfgs), (row - 1) * (row) / 2 + col) << ", ";
                    }
                }
                std::cout << std::endl;
            }
            break;
        } else if (status > 0) {
            std::cout << "Some different error : " << status << std::endl;
            exit(-1);
        } else if (status == -1) {
            std::cout << "Minimization failed, status : " << status << std::endl;
            exit(-1);
        }

        /*for (int row = 0; row < 2016; row++) {
            if (gsl_vector_get(gsl_multimin_fdfminimizer_x(quasi_newton_bfgs), row) < 0) {
                nondiag_neg = true;
                break;
            }
        }
        if (nondiag_neg == true) {

        }*/

        status = gsl_multimin_test_gradient(quasi_newton_bfgs->gradient, std::numeric_limits<double>::epsilon());
        if (status == GSL_SUCCESS) {
            std::cout << "Minimum found in one iteration" << std::endl;
        }
    } while (status == GSL_CONTINUE && iter < 200);

    if (status == GSL_ENOPROG) {
        std::cout << "Skip for a time" << std::endl;
        quasi_newton_bfgs->f = gsl_multimin_fdfminimizer_minimum(quasi_newton_bfgs);
        gsl_vector_memcpy(x, gsl_multimin_fdfminimizer_x(quasi_newton_bfgs));
    }

    double func_value = quasi_newton_bfgs->f;
    gsl_multimin_fdfminimizer_free(quasi_newton_bfgs);
    std::cout << "end maximization step" << std::endl;
    return func_value;
}