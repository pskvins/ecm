//
// Created by 박석환 on 2021/09/02.
//

#include <algorithm>
#include <gsl/gsl_linalg.h>
#include <omp.h>
#include "newick.h"
#include "maximization_step.h"
#include "expectation_step.cpp"

void calculate_derivative(gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, newick_start *start, gsl_matrix *gradient) {
    gsl_matrix_set_all(gradient, 0.0);
    gsl_matrix *dxepon_matrix[64 * 63 / 2];
    for (int num = 0; num < 2016; num++) {
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
        thread_local gsl_matrix *F = gsl_matrix_alloc(64, 64);//freed
        thread_local size_t row_tar = num_to_coordinate.find(diff)->second.begin()->first;
        thread_local size_t col_tar = num_to_coordinate.find(diff)->second.begin()->second;
        while (next_iteration.size() != 1) {
            size_t size = next_iteration.size();
            for (size_t num = 0; num < size; num++) {
                gsl_matrix_set_all(dxepon_matrix[diff], 0.0);
                for (size_t row = 0; row < 64; row++) {
                    for (size_t col = 0; col < 64; col++) {
                        if (row == col) {
                            gsl_matrix_set(dxepon_matrix[diff], row, col,
                                           (2 * codon_freq[row_tar] * codon_freq[col_tar] * diagonal[row] - codon_freq[col_tar] * normalize) / pow(normalize, 2));
                        } else if ((row == row_tar && col == col_tar) || (row == col_tar && col == row_tar)) {
                            gsl_matrix_set(dxepon_matrix[diff], row, col,
                                           (codon_freq[col] * normalize - 2 * codon_freq[row] * pow(codon_freq[col], 2) * gsl_matrix_get(qmatrix, row_tar, col_tar) / pow(normalize, 2)));
                        } else {
                            gsl_matrix_set(dxepon_matrix[diff], row, col,
                                           -gsl_matrix_get(qmatrix, row, col) * 2 * codon_freq[row_tar] * codon_freq[col_tar] / pow(normalize, 2));
                        }
                    }
                }
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigen_inverse, dxepon_matrix[diff], 0.0, dxepon_matrix[diff]);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dxepon_matrix[diff], eigenvector, 0.0, dxepon_matrix[diff]);
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
                gsl_matrix_mul_elements(dxepon_matrix[diff], F);
                gsl_matrix_free(F);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigenvector, dxepon_matrix[diff], 0.0, dxepon_matrix[diff]);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dxepon_matrix[diff], eigen_inverse, 0.0, dxepon_matrix[num]);
                //putting negative of gradient
                for (size_t row = 0; row < 64; row++) {
                    for (size_t col = 0; col < 64; col++) {
                        gsl_matrix_set(gradient, diff, 0, gsl_matrix_get(gradient, diff, 0) - next_iteration[num]->expectation[64 * row + col]
                                                                                                  * gsl_matrix_get(dxepon_matrix[diff], row, col)
                                                                                                  / gsl_matrix_get(next_iteration[num]->expon_matrix[0], row, col));
                    }
                }
                if (std::find(next_iteration.begin(), next_iteration.end(), next_iteration[num]->next) == next_iteration.end()) {
                    if (next_iteration[num]->next != NULL) {
                        next_iteration.emplace_back(next_iteration[num]->next);
                    }
                }
            }
            next_iteration.erase(next_iteration.begin(), next_iteration.begin() + size);
        }
    }
    gsl_matrix_free(eigenvector);
    gsl_matrix_free(eigen_inverse);
    gsl_vector_free(eigenvalue);
    for (int num = 0; num < 2016; num++) {
        gsl_matrix_free(dxepon_matrix[num]);
    }
}

void normalize_qmatirx(gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *direction, double lambda, gsl_matrix *new_matrix) {
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            gsl_matrix_set(new_matrix, row, col, (gsl_matrix_get(qmatrix, row, col) + lambda * gsl_matrix_get(direction, 64 * row + col, 0)) * codon_freq[col]);
        }
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
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            gsl_matrix_set(new_matrix, row, col, gsl_matrix_get(new_matrix, row, col)/sum);
        }
    }
}

double calculate_maximization_function(newick_start *start, short int index) {
    double result = 0;
    std::vector<newick_graph*> next_iteration = start->next;
    while(next_iteration.size() != 1) {
        size_t size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            for (size_t row = 0; row < 64; row++) {
                for (size_t col = 0; col < 64; col++) {
                    result -= next_iteration[num]->expectation[64 * row + col] * log(gsl_matrix_get(next_iteration[num]->expon_matrix[index], row, col));
                }
            }
            if (std::find(next_iteration.begin(), next_iteration.end(), next_iteration[num]->next) == next_iteration.end()) {
                if (next_iteration[num]->next != NULL) {
                    next_iteration.emplace_back(next_iteration[num]->next);
                }
            }
        }
        next_iteration.erase(next_iteration.begin(), next_iteration.begin() + size);
    }
    return result;
}

void update_expon_matrix(newick_start *start, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue) {
    std::vector<newick_graph*> next_iterator = start->next;
    gsl_matrix *eigenvector_temp = gsl_matrix_alloc(64, 64);//freed
    gsl_vector *eigenvalue_temp = gsl_vector_alloc(64);//freed
    while (next_iterator.size() != 1) {
        size_t size = next_iterator.size();
        for (size_t num = 0; num < size; num++) {
            gsl_matrix_memcpy(next_iterator[num]->expon_matrix[1], calculate_expon_matrix(eigenvector, eigen_inverse, eigenvalue, eigenvector_temp, eigenvalue_temp, next_iterator[num]->branch_length));
            if (std::find(next_iterator.begin(), next_iterator.end(), next_iterator[num]->next) == next_iterator.end()) {
                if (next_iterator[num]->next != NULL) {
                    next_iterator.emplace_back(next_iterator[num]->next);
                }
            }
        }
        next_iterator.erase(next_iterator.begin(), next_iterator.begin() + size);
    }
    gsl_matrix_free(eigenvector_temp);
    gsl_vector_free(eigenvalue_temp);
}

double backtracking(newick_start *start, gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *direction, gsl_matrix *gradient, double &lamb, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue) {
    gsl_matrix *qmatrix_temp = gsl_matrix_alloc(64, 64);//freed
    gsl_matrix_memcpy(qmatrix_temp, qmatrix);
    //TODO: Do I have to normalize in this step??-->seems like has to be normalized (equation needs normalization)
    double funcg_ori = calculate_maximization_function(start, 0);
    double slope = 0.0;
    for (int num = 0; num < 2016; num++) {
        slope += gsl_matrix_get(gradient, num, 0) * gsl_matrix_get(direction, num, 0);
    }
    assert(slope < 0.0 && "Direction not going down");
    const double machine_precision = 2.0e-14;
    const double alpha = 1.0e-4;
    double temp;
    double test = 0.0;
    size_t row;
    size_t col;
    for (size_t num = 0; num < 63 * 64 / 2; num++) {
        row = num_to_coordinate.find(num)->second.begin()->first;
        col = num_to_coordinate.find(num)->second.begin()->second;
        temp = abs(gsl_matrix_get(direction, num, 0)) / std::max(abs(gsl_matrix_get(qmatrix, row, col)), 1.0);
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
        //normalize after adding direction
        normalize_qmatirx(qmatrix, codon_freq, direction, lamb, qmatrix_temp);

        //get eigenvalues and eigenvectors
        gsl_matrix *qmatrix_cal = gsl_matrix_alloc(64, 64);//freed
        get_eigenvector_and_inverse(qmatrix_temp, codon_freq, qmatrix_cal, eigenvector, eigen_inverse, eigenvalue);

        //update matrix in tree
        update_expon_matrix(start, eigenvector, eigen_inverse, eigenvalue);

        //calculate g(lambda)
        funcg = calculate_maximization_function(start, 1);
        if (lamb < lambmin) {
            gsl_matrix_free(qmatrix_temp);
            lamb = 0.0;
            return funcg;
        } else if (funcg <= funcg_ori + alpha * lamb * slope) {
            gsl_matrix_free(qmatrix_temp);
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
                    if (inside_sqrt < 0.0) {
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

double quasi_Newton_method(newick_start *start, gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue) {
    const double machine_precision = 2.0e-14;
    const int iterate_max = 200;
    //initialize matrix
    gsl_matrix *hessian = gsl_matrix_alloc(2016, 2016);//freed
    gsl_matrix *hessian_append = gsl_matrix_alloc(2016, 2016);//freed
    gsl_matrix_set_all(hessian, 0.0);
    gsl_matrix_set_all(hessian_append, 0.0);
    for (int num = 0; num < 2016; num++) {
        gsl_matrix_set(hessian, num, num, 1.0);
    }
    gsl_matrix *gradient = gsl_matrix_alloc(2016, 1);//freed
    gsl_matrix *gradient_new = gsl_matrix_alloc(2016, 1);//freed
    gsl_matrix_set_all(gradient, 0.0);
    gsl_matrix_set_all(gradient_new, 0.0);
    calculate_derivative(qmatrix, codon_freq, eigenvector, eigen_inverse, eigenvalue, start, gradient);
    gsl_matrix *dx = gsl_matrix_alloc(2016, 1);//freed
    gsl_matrix *x = gsl_matrix_alloc(2016, 1);//freed
    gsl_matrix *x_new = gsl_matrix_alloc(2016, 1);//freed
    for (int num = 0; num < 2016; num++) {
        gsl_matrix_set(dx, num, 0, -gsl_matrix_get(gradient, num, 0));
    }
    int row_map, col_map;
    for (int num = 0; num < 2016; num++) {
        row_map = num_to_coordinate.find(num)->second.begin()->first;
        col_map = num_to_coordinate.find(num)->second.begin()->second;
        gsl_matrix_set(x, num, 0, gsl_matrix_get(qmatrix, row_map, col_map));
    }

    gsl_matrix *qmatrix_new = gsl_matrix_alloc(64, 64);//freed
    double lambda = 1.0;
    double function_value;
    //start loop
    for (int its = 0; its < iterate_max; its++) {
        function_value = backtracking(start, qmatrix, codon_freq, dx, gradient, lambda, eigenvector, eigen_inverse, eigenvalue);
        for (int num = 0; num < 2016; num++) {
            gsl_matrix_set(dx, num, 0, lambda * gsl_matrix_get(dx, num, 0));
            row_map = num_to_coordinate.find(num)->second.begin()->first;
            col_map = num_to_coordinate.find(num)->second.begin()->second;
            gsl_matrix_set(qmatrix_new, row_map, col_map, gsl_matrix_get(qmatrix, row_map, col_map) + gsl_matrix_get(dx, num, 0));
            gsl_matrix_set(qmatrix_new, col_map, row_map, gsl_matrix_get(qmatrix_new, row_map, col_map));
            gsl_matrix_set(x_new, num, 0, gsl_matrix_get(x, num, 0) + gsl_matrix_get(dx, num, 0));
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
            double temp = abs(gsl_matrix_get(dx, num, 0)) / std::max(abs(gsl_matrix_get(x_new, num, 0)), 1.0);
            if (temp > test) {
                test = temp;
            }
        }
        if (test < 4 * machine_precision) {
            gsl_matrix_free(hessian);
            gsl_matrix_free(hessian_append);
            gsl_matrix_free(gradient);
            gsl_matrix_free(gradient_new);
            gsl_matrix_free(x);
            gsl_matrix_free(x_new);
            gsl_matrix_free(dx);
            gsl_matrix_free(qmatrix_new);
            return function_value;
        }
        calculate_derivative(qmatrix_new, codon_freq, eigenvector, eigen_inverse, eigenvalue, start, gradient_new);
        //now gradient is the yk
        for (int num = 0 ; num < 2016; num++) {
            gsl_matrix_set(gradient, num, 0, gsl_matrix_get(gradient_new, num, 0) - gsl_matrix_get(gradient, num, 0));
        }
        //update hessian matrix
        double denominator = 0;
        for (int num = 0; num < 2016; num++) {
            denominator += gsl_matrix_get(gradient, num, 0) * gsl_matrix_get(x_new, num, 0);
        }
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0 / denominator, dx, gradient, 0.0, hessian_append);
        for (int num = 0; num < 2016; num++) {
            gsl_matrix_set(hessian_append, num, num, 1.0 + gsl_matrix_get(hessian_append, num, num));
        }
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, hessian_append, hessian, 0.0, hessian);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0 / denominator, gradient, dx, 0.0, hessian_append);
        for (int num = 0; num < 2016; num++) {
            gsl_matrix_set(hessian_append, num, num, 1.0 + gsl_matrix_get(hessian_append, num, num));
        }
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, hessian, hessian_append, 0.0, hessian);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / denominator, dx, dx, 0.0, hessian_append);
        gsl_matrix_add(hessian, hessian_append);
        //update gradient, qmatirx x and dx
        gsl_matrix_memcpy(x, x_new);
        gsl_matrix_memcpy(gradient, gradient_new);
        gsl_matrix_memcpy(qmatrix, qmatrix_new);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, hessian, gradient, 0.0, dx);
    }
    throw("Too many iterations");
}