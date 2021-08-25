//
// Created by 박석환 on 2021/08/13.
//

#include <vector>
#include <algorithm>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>


double cal_prob_from_qmatrix(gsl_matrix *qmatrix, double codon_freq, char base_before, char base_after, float branch_length) {
    //get eigenvalues and left & right eigenvectors
    //set memories
    gsl_matrix *original_matrix = gsl_matrix_alloc(64, 64);
    gsl_matrix_memcpy(original_matrix, qmatrix);
    gsl_matrix *qmatrix_temp = gsl_matrix_alloc(64, 64);//freed
    gsl_matrix_memcpy(qmatrix_temp, qmatrix);
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            assert(gsl_matrix_get(original_matrix, row, col) == gsl_matrix_get(qmatrix_temp, row, col));
        }
    }
    gsl_vector_complex *eigenvalue = gsl_vector_complex_alloc(64);
    gsl_vector_complex *eigenvalue_t = gsl_vector_complex_alloc(64);//freed
    gsl_matrix_complex *right_eigen = gsl_matrix_complex_alloc(64, 64); //freed
    gsl_matrix_complex *left_eigen = gsl_matrix_complex_alloc(64, 64); //freed
    gsl_eigen_nonsymmv_workspace *eigen_space = gsl_eigen_nonsymmv_alloc(64);//freed

    // calculate right eigenvector
    gsl_eigen_nonsymmv(qmatrix_temp, eigenvalue, right_eigen, eigen_space);

    //now do the same with transposed qmatrix
    gsl_matrix_transpose_memcpy(qmatrix_temp, qmatrix);
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            assert(gsl_matrix_get(original_matrix, row, col) == gsl_matrix_get(qmatrix_temp, col, row));
        }
    }

    // calculate left eigenvector (transposed)
    gsl_eigen_nonsymmv(qmatrix_temp, eigenvalue_t, left_eigen, eigen_space);
    gsl_matrix_free(qmatrix_temp);
    gsl_eigen_nonsymmv_free(eigen_space);

    //sort right & left eigenvector so that they can have same order
    gsl_eigen_nonsymmv_sort(eigenvalue, right_eigen, GSL_EIGEN_SORT_ABS_ASC);
    gsl_eigen_nonsymmv_sort(eigenvalue_t, left_eigen, GSL_EIGEN_SORT_ABS_ASC);

    //fix mechanical error (make 0 if the value is smaller than 1.0e-15)
    gsl_complex check;
    gsl_complex check_t;
    gsl_complex set_data;
    for (size_t i = 0; i < 64; i++) {
        check = gsl_vector_complex_get(eigenvalue, i);
        check_t = gsl_vector_complex_get(eigenvalue_t, i);
        if (abs(check.dat[0]) < 1.0e-15 && check.dat[0] != 0) {
            set_data = gsl_complex_rect(0.0,check.dat[1]);
            gsl_vector_complex_set(eigenvalue, i, set_data);
        }
        if (abs(check.dat[1]) < 1.0e-15 && check.dat[1] != 0) {
            set_data = gsl_complex_rect(check.dat[0], 0.0);
            gsl_vector_complex_set(eigenvalue, i, set_data);
        }
        if (abs(check_t.dat[0]) < 1.0e-15 && check_t.dat[0] != 0) {
            set_data = gsl_complex_rect(0.0, check_t.dat[1]);
            gsl_vector_complex_set(eigenvalue_t, i, set_data);
        }
        if (abs(check_t.dat[1]) < 1.0e-15 && check_t.dat[1] != 0) {
            set_data = gsl_complex_rect(check_t.dat[0], 0.0);
            gsl_vector_complex_set(eigenvalue_t, i, set_data);
        }
    }

    //check if all the eigenvalues are real number && eigenvalue == eigenvalue_t
    for (size_t i = 0; i < eigenvalue->size; i++) {
        check = gsl_vector_complex_get(eigenvalue, i);
        check_t = gsl_vector_complex_get(eigenvalue_t, i);
        assert(abs(check.dat[0] - check_t.dat[0]) < 1.0e-13);
        assert(abs(check.dat[1] - check_t.dat[1]) < 1.0e-13);
        assert(check.dat[1] == 0);
    }


    std::cout << "right eigenvector" << std::endl;
    //for debug, check eigenvectors
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            check = gsl_matrix_complex_get(right_eigen, col, row);
            if (check.dat[1] != 0) {
                std::cout << check.dat[0] << " " << check.dat[1] << std::endl;
            }
        }
        std::cout << "new column" << std::endl;
    }

    std::cout << "left eigenvector" << std::endl;
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            check_t= gsl_matrix_complex_get(left_eigen, col, row);
            if (check_t.dat[1] != 0) {
            std::cout << check_t.dat[0] << " " << check_t.dat[1] << std::endl;
            }
        }
        std::cout << "new row" << std::endl;
    }

    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            check = gsl_matrix_complex_get(right_eigen, row, col);
            check_t = gsl_matrix_complex_get(left_eigen, row, col);
            assert(check.dat[1] == 0 && check_t.dat[1] == 0);
        }
    }

    gsl_vector_complex_free(eigenvalue_t);
    //make left eigenvector transversed
    gsl_matrix_complex *left_eigen_t = gsl_matrix_complex_alloc(64, 64); //freed
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            gsl_matrix_complex_set(left_eigen_t, row, col, gsl_matrix_complex_get(left_eigen, col, row));
        }
    }
    gsl_matrix_complex_free(left_eigen);
    //check if all eigenvalues are real number

    gsl_matrix *right_eigen_real = gsl_matrix_alloc(64, 64);
    gsl_matrix *left_eigen_real = gsl_matrix_alloc(64, 64);
    for (size_t row = 0; row < 64; row++) {
        for (size_t col = 0; col < 64; col++) {
            gsl_matrix_set(right_eigen_real, row, col, gsl_matrix_complex_get(right_eigen, row, col).dat[0]);
            gsl_matrix_set(left_eigen_real, row, col, gsl_matrix_complex_get(left_eigen_t, row, col).dat[0]);
        }
    }
    gsl_matrix_complex_free(right_eigen);
    gsl_matrix_complex_free(left_eigen_t);

    gsl_vector *expon_eigenvalue = gsl_vector_alloc(64);

    //calculate possibility
    for (size_t i = 0; i < 64; i++) {
        if (gsl_vector_complex_get(eigenvalue, i).dat[0] != 0) {
            gsl_vector_set(expon_eigenvalue, i, pow(M_E, gsl_vector_complex_get(eigenvalue, i).dat[0]) * branch_length);
        }
        else {
            gsl_vector_set(expon_eigenvalue, i, 0.0);
        }
    }

    return 0.0;
}