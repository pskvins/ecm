//
// Created by 박석환 on 2021/08/13.
//

#include <algorithm>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "process_maf.h"

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

gsl_matrix* calculate_expon_matrix(gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue,
                                  gsl_matrix *eigenvector_temp, gsl_vector *eigenvalue_temp, double branch_length) {
    for (size_t num = 0; num < 64; num++) {
        gsl_vector_set(eigenvalue_temp, num, pow(M_E, gsl_vector_get(eigenvalue, num) * branch_length));
    }
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            gsl_matrix_set(eigenvector_temp, row, col, gsl_matrix_get(eigenvector, row, col) * gsl_vector_get(eigenvalue_temp, col));
        }
    }
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigenvector_temp, eigen_inverse, 0.0, eigenvector_temp);
    return eigenvector_temp;
}

void set_matrices(newick_start *start, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue) {
    std::vector<newick_graph*> next_iterator = start->next;
    gsl_matrix *eigenvector_temp = gsl_matrix_alloc(64, 64);//freed
    gsl_vector *eigenvalue_temp = gsl_vector_alloc(64);//freed
    while (next_iterator.size() != 1) {
        size_t size = next_iterator.size();
        for (size_t num = 0; num < size; num++) {
            gsl_matrix_memcpy(next_iterator[num]->expon_matrix,calculate_expon_matrix(eigenvector, eigen_inverse, eigenvalue, eigenvector_temp, eigenvalue_temp, next_iterator[num]->branch_length));
        }
    }
    gsl_matrix_free(eigenvector_temp);
    gsl_vector_free(eigenvalue_temp);
}

double felsenstein_algorithm(newick_graph *node, char base, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue) {
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
            sum_first += node->previous[0]->felsenstein[num] * gsl_matrix_get(node->previous[0]->expon_matrix, base, num);
            sum_second += node->previous[1]->felsenstein[num] * gsl_matrix_get(node->previous[1]->expon_matrix, base, num);
        }
        node->felsenstein[base] = sum_first * sum_second;
        return sum_first * sum_second;
    }
}

void conduct_felsenstein(newick_start *start, aligned_codon *codon_set, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue) {
    std::vector<newick_graph*> next_iteration = start->next;
    for (size_t num = 0; num < next_iteration.size(); num++) {
        for (size_t set = 0; set < sizeof(codon_set) / sizeof(codon_set[0]); set++) {
            if (next_iteration[num]->species == codon_set[set].species) {
                for (size_t i = 0; i < 64; i++) {
                    if (i != codon_set[set].codon) {
                        next_iteration[num]->base[i] = false;
                    }
                }
                break;
            }
        }
    }
    while (!next_iteration.empty()) {
        int size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            for (char codon = 0; codon < 64; codon++) {
                felsenstein_algorithm(next_iteration[num], codon, eigenvector, eigen_inverse, eigenvalue);
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

void calculate_upper(newick_graph *end, double *codon_freq) {
    std::vector<newick_graph*> next_iteration = end->previous;
    for (size_t num = 0; num < 64; num++) {
        end->upper[num] = codon_freq[num];
    }
    while (!next_iteration.empty()) {
        int size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            for (char base = 0; base < 64; base++) {
                double fel_upper = 0;
                for (char base_first = 0; base_first < 64; base_first++) {
                    for (char base_second = 0; base_second < 64; base_second++) {
                        fel_upper += next_iteration[num]->next->upper[base_first] *
                                     gsl_matrix_get(next_iteration[num]->expon_matrix, base_first, base)
                                     * next_iteration[num]->next->previous[(next_iteration[num] == next_iteration[num]->next->previous[0]) ? 1 : 0]->felsenstein[base_second]
                                     * gsl_matrix_get(next_iteration[num]->next->previous[(next_iteration[num] == next_iteration[num]->next->previous[0]) ? 1 : 0]->expon_matrix, base_first, base_second);
                        if (!next_iteration[num]->previous.empty()) {
                            next_iteration.emplace_back(next_iteration[num]->previous[0]);
                            next_iteration.emplace_back(next_iteration[num]->previous[1]);
                        }
                    }
                }
            }
        }
        next_iteration.erase(next_iteration.begin(), next_iteration.begin() + size);
    }
}

void update_upper(newick_start *start, double *codon_freq) {
    //calculate denominator
    double denominator = 1.0;
    std::vector<newick_graph*> next_iteration = start->next;
    for (size_t num = 0; num < next_iteration.size(); num++) {
        for (char base = 0; base < 64; base++) {
            if (next_iteration[num]->base[base] == true) {
                denominator *= codon_freq[base];
                break;
            }
        }
    }

    //update upper
    while (!next_iteration.empty()) {
        int size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            for (char base = 0; base < 64; base++) {
                next_iteration[num]->upper[base] = next_iteration[num]->felsenstein[base] * next_iteration[num]->upper[base] / denominator;
            }
            if (next_iteration[num]->next != NULL) {
                next_iteration.emplace_back(next_iteration[num]->next);
            }
        }
        next_iteration.erase(next_iteration.begin(), next_iteration.begin() + size);
    }
}

void calculate_expectation(newick_start *start) {
    std::vector<newick_graph*> next_iteration = start->next;
    while (next_iteration.size() != 1) {
        int size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            for (char base_next = 0; base_next < 64; base_next++) {
                double denominator = 0.0;
                for (char base_acc = 0; base_acc < 64; base_acc++) {
                    denominator += next_iteration[num]->felsenstein[base_acc] * gsl_matrix_get(next_iteration[num]->expon_matrix, base_next, base_acc);
                }
                for (char base_curr = 0; base_curr < 64; base_curr++) {
                    next_iteration[num]->expectation[64 * base_next + base_curr] += next_iteration[num]->felsenstein[base_curr] *
                            gsl_matrix_get(next_iteration[num]->expon_matrix, base_next, base_curr) * next_iteration[num]->next->upper[base_next] / denominator;
                }
            }
            if (next_iteration[num]->next != NULL) {
                next_iteration.emplace_back(next_iteration[num]->next);
            }
        }
        next_iteration.erase(next_iteration.begin(), next_iteration.begin() + size);
    }
}

//Integrating all the expectation steps
void conduct_expectation_step(std::vector<aligned_codon> aligned_codon_set, newick_start *start, newick_graph *end,
                              gsl_matrix *qmatrix, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, double *codon_freq) {
    std::vector<newick_graph*> next_iterator = start->next;
    for (size_t num = 0; num < next_iterator.size(); num++) {
        for (char base = 0; base < 64; base++) {
            next_iterator[num]->base[base] = true;
        }
    }

    while (!next_iterator.empty()) {
        size_t size = next_iterator.size();
        for (size_t num = 0; num < size; num++) {
            for (char base_first = 0; base_first < 64; base_first++) {
                for (char base_second = 0; base_second < 64; base_second++) {
                    next_iterator[num]->expectation[64 * base_first + base_second] = 0.0;
                }
            }
            if (next_iterator[num]->next != NULL) {
                next_iterator.emplace_back(next_iterator[num]->next);
            }
        }
        next_iterator.erase(next_iterator.begin(), next_iterator.begin() + size);
    }

    get_eigenvector_and_inverse(qmatrix, 0, eigenvector, eigen_inverse, eigenvalue);

    set_matrices(start, eigenvector, eigen_inverse, eigenvalue);

    for (size_t num = 0; num < aligned_codon_set.size(); num++) {
        conduct_felsenstein(start, &aligned_codon_set[num], eigenvector, eigen_inverse, eigenvalue);

        calculate_upper(end, codon_freq);

        update_upper(start, codon_freq);

        calculate_expectation(start);
    }
}