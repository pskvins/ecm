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
#include <assert.h>
#include <omp.h>

//map for maximization step
double num_to_coordinate[2016][2] = {{1, 0},
                                     {2, 0},{2, 1},
                                     {3, 0},{3, 1},{3, 2},
                                     {4, 0},{4, 1},{4, 2},{4, 3},
                                     {5, 0},{5, 1},{5, 2},{5, 3},{5, 4},
                                     {6, 0},{6, 1},{6, 2},{6, 3},{6, 4},{6, 5},
                                     {7, 0},{7, 1},{7, 2},{7, 3},{7, 4},{7, 5},{7, 6},
                                     {8, 0},{8, 1},{8, 2},{8, 3},{8, 4},{8, 5},{8, 6},{8, 7},
                                     {9, 0},{9, 1},{9, 2},{9, 3},{9, 4},{9, 5},{9, 6},{9, 7},{9, 8},
                                     {10, 0},{10, 1},{10, 2},{10, 3},{10, 4},{10, 5},{10, 6},{10, 7},{10, 8},{10, 9},
                                     {11, 0},{11, 1},{11, 2},{11, 3},{11, 4},{11, 5},{11, 6},{11, 7},{11, 8},{11, 9},{11, 10},
                                     {12, 0},{12, 1},{12, 2},{12, 3},{12, 4},{12, 5},{12, 6},{12, 7},{12, 8},{12, 9},{12, 10},{12, 11},
                                     {13, 0},{13, 1},{13, 2},{13, 3},{13, 4},{13, 5},{13, 6},{13, 7},{13, 8},{13, 9},{13, 10},{13, 11},{13, 12},
                                     {14, 0},{14, 1},{14, 2},{14, 3},{14, 4},{14, 5},{14, 6},{14, 7},{14, 8},{14, 9},{14, 10},{14, 11},{14, 12},{14, 13},
                                     {15, 0},{15, 1},{15, 2},{15, 3},{15, 4},{15, 5},{15, 6},{15, 7},{15, 8},{15, 9},{15, 10},{15, 11},{15, 12},{15, 13},{15, 14},
                                     {16, 0},{16, 1},{16, 2},{16, 3},{16, 4},{16, 5},{16, 6},{16, 7},{16, 8},{16, 9},{16, 10},{16, 11},{16, 12},{16, 13},{16, 14},{16, 15},
                                     {17, 0},{17, 1},{17, 2},{17, 3},{17, 4},{17, 5},{17, 6},{17, 7},{17, 8},{17, 9},{17, 10},{17, 11},{17, 12},{17, 13},{17, 14},{17, 15},{17, 16},
                                     {18, 0},{18, 1},{18, 2},{18, 3},{18, 4},{18, 5},{18, 6},{18, 7},{18, 8},{18, 9},{18, 10},{18, 11},{18, 12},{18, 13},{18, 14},{18, 15},{18, 16},{18, 17},
                                     {19, 0},{19, 1},{19, 2},{19, 3},{19, 4},{19, 5},{19, 6},{19, 7},{19, 8},{19, 9},{19, 10},{19, 11},{19, 12},{19, 13},{19, 14},{19, 15},{19, 16},{19, 17},{19, 18},
                                     {20, 0},{20, 1},{20, 2},{20, 3},{20, 4},{20, 5},{20, 6},{20, 7},{20, 8},{20, 9},{20, 10},{20, 11},{20, 12},{20, 13},{20, 14},{20, 15},{20, 16},{20, 17},{20, 18},{20, 19},
                                     {21, 0},{21, 1},{21, 2},{21, 3},{21, 4},{21, 5},{21, 6},{21, 7},{21, 8},{21, 9},{21, 10},{21, 11},{21, 12},{21, 13},{21, 14},{21, 15},{21, 16},{21, 17},{21, 18},{21, 19},{21, 20},
                                     {22, 0},{22, 1},{22, 2},{22, 3},{22, 4},{22, 5},{22, 6},{22, 7},{22, 8},{22, 9},{22, 10},{22, 11},{22, 12},{22, 13},{22, 14},{22, 15},{22, 16},{22, 17},{22, 18},{22, 19},{22, 20},{22, 21},
                                     {23, 0},{23, 1},{23, 2},{23, 3},{23, 4},{23, 5},{23, 6},{23, 7},{23, 8},{23, 9},{23, 10},{23, 11},{23, 12},{23, 13},{23, 14},{23, 15},{23, 16},{23, 17},{23, 18},{23, 19},{23, 20},{23, 21},{23, 22},
                                     {24, 0},{24, 1},{24, 2},{24, 3},{24, 4},{24, 5},{24, 6},{24, 7},{24, 8},{24, 9},{24, 10},{24, 11},{24, 12},{24, 13},{24, 14},{24, 15},{24, 16},{24, 17},{24, 18},{24, 19},{24, 20},{24, 21},{24, 22},{24, 23},
                                     {25, 0},{25, 1},{25, 2},{25, 3},{25, 4},{25, 5},{25, 6},{25, 7},{25, 8},{25, 9},{25, 10},{25, 11},{25, 12},{25, 13},{25, 14},{25, 15},{25, 16},{25, 17},{25, 18},{25, 19},{25, 20},{25, 21},{25, 22},{25, 23},{25, 24},
                                     {26, 0},{26, 1},{26, 2},{26, 3},{26, 4},{26, 5},{26, 6},{26, 7},{26, 8},{26, 9},{26, 10},{26, 11},{26, 12},{26, 13},{26, 14},{26, 15},{26, 16},{26, 17},{26, 18},{26, 19},{26, 20},{26, 21},{26, 22},{26, 23},{26, 24},{26, 25},
                                     {27, 0},{27, 1},{27, 2},{27, 3},{27, 4},{27, 5},{27, 6},{27, 7},{27, 8},{27, 9},{27, 10},{27, 11},{27, 12},{27, 13},{27, 14},{27, 15},{27, 16},{27, 17},{27, 18},{27, 19},{27, 20},{27, 21},{27, 22},{27, 23},{27, 24},{27, 25},{27, 26},
                                     {28, 0},{28, 1},{28, 2},{28, 3},{28, 4},{28, 5},{28, 6},{28, 7},{28, 8},{28, 9},{28, 10},{28, 11},{28, 12},{28, 13},{28, 14},{28, 15},{28, 16},{28, 17},{28, 18},{28, 19},{28, 20},{28, 21},{28, 22},{28, 23},{28, 24},{28, 25},{28, 26},{28, 27},
                                     {29, 0},{29, 1},{29, 2},{29, 3},{29, 4},{29, 5},{29, 6},{29, 7},{29, 8},{29, 9},{29, 10},{29, 11},{29, 12},{29, 13},{29, 14},{29, 15},{29, 16},{29, 17},{29, 18},{29, 19},{29, 20},{29, 21},{29, 22},{29, 23},{29, 24},{29, 25},{29, 26},{29, 27},{29, 28},
                                     {30, 0},{30, 1},{30, 2},{30, 3},{30, 4},{30, 5},{30, 6},{30, 7},{30, 8},{30, 9},{30, 10},{30, 11},{30, 12},{30, 13},{30, 14},{30, 15},{30, 16},{30, 17},{30, 18},{30, 19},{30, 20},{30, 21},{30, 22},{30, 23},{30, 24},{30, 25},{30, 26},{30, 27},{30, 28},{30, 29},
                                     {31, 0},{31, 1},{31, 2},{31, 3},{31, 4},{31, 5},{31, 6},{31, 7},{31, 8},{31, 9},{31, 10},{31, 11},{31, 12},{31, 13},{31, 14},{31, 15},{31, 16},{31, 17},{31, 18},{31, 19},{31, 20},{31, 21},{31, 22},{31, 23},{31, 24},{31, 25},{31, 26},{31, 27},{31, 28},{31, 29},{31, 30},
                                     {32, 0},{32, 1},{32, 2},{32, 3},{32, 4},{32, 5},{32, 6},{32, 7},{32, 8},{32, 9},{32, 10},{32, 11},{32, 12},{32, 13},{32, 14},{32, 15},{32, 16},{32, 17},{32, 18},{32, 19},{32, 20},{32, 21},{32, 22},{32, 23},{32, 24},{32, 25},{32, 26},{32, 27},{32, 28},{32, 29},{32, 30},{32, 31},
                                     {33, 0},{33, 1},{33, 2},{33, 3},{33, 4},{33, 5},{33, 6},{33, 7},{33, 8},{33, 9},{33, 10},{33, 11},{33, 12},{33, 13},{33, 14},{33, 15},{33, 16},{33, 17},{33, 18},{33, 19},{33, 20},{33, 21},{33, 22},{33, 23},{33, 24},{33, 25},{33, 26},{33, 27},{33, 28},{33, 29},{33, 30},{33, 31},{33, 32},
                                     {34, 0},{34, 1},{34, 2},{34, 3},{34, 4},{34, 5},{34, 6},{34, 7},{34, 8},{34, 9},{34, 10},{34, 11},{34, 12},{34, 13},{34, 14},{34, 15},{34, 16},{34, 17},{34, 18},{34, 19},{34, 20},{34, 21},{34, 22},{34, 23},{34, 24},{34, 25},{34, 26},{34, 27},{34, 28},{34, 29},{34, 30},{34, 31},{34, 32},{34, 33},
                                     {35, 0},{35, 1},{35, 2},{35, 3},{35, 4},{35, 5},{35, 6},{35, 7},{35, 8},{35, 9},{35, 10},{35, 11},{35, 12},{35, 13},{35, 14},{35, 15},{35, 16},{35, 17},{35, 18},{35, 19},{35, 20},{35, 21},{35, 22},{35, 23},{35, 24},{35, 25},{35, 26},{35, 27},{35, 28},{35, 29},{35, 30},{35, 31},{35, 32},{35, 33},{35, 34},
                                     {36, 0},{36, 1},{36, 2},{36, 3},{36, 4},{36, 5},{36, 6},{36, 7},{36, 8},{36, 9},{36, 10},{36, 11},{36, 12},{36, 13},{36, 14},{36, 15},{36, 16},{36, 17},{36, 18},{36, 19},{36, 20},{36, 21},{36, 22},{36, 23},{36, 24},{36, 25},{36, 26},{36, 27},{36, 28},{36, 29},{36, 30},{36, 31},{36, 32},{36, 33},{36, 34},{36, 35},
                                     {37, 0},{37, 1},{37, 2},{37, 3},{37, 4},{37, 5},{37, 6},{37, 7},{37, 8},{37, 9},{37, 10},{37, 11},{37, 12},{37, 13},{37, 14},{37, 15},{37, 16},{37, 17},{37, 18},{37, 19},{37, 20},{37, 21},{37, 22},{37, 23},{37, 24},{37, 25},{37, 26},{37, 27},{37, 28},{37, 29},{37, 30},{37, 31},{37, 32},{37, 33},{37, 34},{37, 35},{37, 36},
                                     {38, 0},{38, 1},{38, 2},{38, 3},{38, 4},{38, 5},{38, 6},{38, 7},{38, 8},{38, 9},{38, 10},{38, 11},{38, 12},{38, 13},{38, 14},{38, 15},{38, 16},{38, 17},{38, 18},{38, 19},{38, 20},{38, 21},{38, 22},{38, 23},{38, 24},{38, 25},{38, 26},{38, 27},{38, 28},{38, 29},{38, 30},{38, 31},{38, 32},{38, 33},{38, 34},{38, 35},{38, 36},{38, 37},
                                     {39, 0},{39, 1},{39, 2},{39, 3},{39, 4},{39, 5},{39, 6},{39, 7},{39, 8},{39, 9},{39, 10},{39, 11},{39, 12},{39, 13},{39, 14},{39, 15},{39, 16},{39, 17},{39, 18},{39, 19},{39, 20},{39, 21},{39, 22},{39, 23},{39, 24},{39, 25},{39, 26},{39, 27},{39, 28},{39, 29},{39, 30},{39, 31},{39, 32},{39, 33},{39, 34},{39, 35},{39, 36},{39, 37},{39, 38},
                                     {40, 0},{40, 1},{40, 2},{40, 3},{40, 4},{40, 5},{40, 6},{40, 7},{40, 8},{40, 9},{40, 10},{40, 11},{40, 12},{40, 13},{40, 14},{40, 15},{40, 16},{40, 17},{40, 18},{40, 19},{40, 20},{40, 21},{40, 22},{40, 23},{40, 24},{40, 25},{40, 26},{40, 27},{40, 28},{40, 29},{40, 30},{40, 31},{40, 32},{40, 33},{40, 34},{40, 35},{40, 36},{40, 37},{40, 38},{40, 39},
                                     {41, 0},{41, 1},{41, 2},{41, 3},{41, 4},{41, 5},{41, 6},{41, 7},{41, 8},{41, 9},{41, 10},{41, 11},{41, 12},{41, 13},{41, 14},{41, 15},{41, 16},{41, 17},{41, 18},{41, 19},{41, 20},{41, 21},{41, 22},{41, 23},{41, 24},{41, 25},{41, 26},{41, 27},{41, 28},{41, 29},{41, 30},{41, 31},{41, 32},{41, 33},{41, 34},{41, 35},{41, 36},{41, 37},{41, 38},{41, 39},{41, 40},
                                     {42, 0},{42, 1},{42, 2},{42, 3},{42, 4},{42, 5},{42, 6},{42, 7},{42, 8},{42, 9},{42, 10},{42, 11},{42, 12},{42, 13},{42, 14},{42, 15},{42, 16},{42, 17},{42, 18},{42, 19},{42, 20},{42, 21},{42, 22},{42, 23},{42, 24},{42, 25},{42, 26},{42, 27},{42, 28},{42, 29},{42, 30},{42, 31},{42, 32},{42, 33},{42, 34},{42, 35},{42, 36},{42, 37},{42, 38},{42, 39},{42, 40},{42, 41},
                                     {43, 0},{43, 1},{43, 2},{43, 3},{43, 4},{43, 5},{43, 6},{43, 7},{43, 8},{43, 9},{43, 10},{43, 11},{43, 12},{43, 13},{43, 14},{43, 15},{43, 16},{43, 17},{43, 18},{43, 19},{43, 20},{43, 21},{43, 22},{43, 23},{43, 24},{43, 25},{43, 26},{43, 27},{43, 28},{43, 29},{43, 30},{43, 31},{43, 32},{43, 33},{43, 34},{43, 35},{43, 36},{43, 37},{43, 38},{43, 39},{43, 40},{43, 41},{43, 42},
                                     {44, 0},{44, 1},{44, 2},{44, 3},{44, 4},{44, 5},{44, 6},{44, 7},{44, 8},{44, 9},{44, 10},{44, 11},{44, 12},{44, 13},{44, 14},{44, 15},{44, 16},{44, 17},{44, 18},{44, 19},{44, 20},{44, 21},{44, 22},{44, 23},{44, 24},{44, 25},{44, 26},{44, 27},{44, 28},{44, 29},{44, 30},{44, 31},{44, 32},{44, 33},{44, 34},{44, 35},{44, 36},{44, 37},{44, 38},{44, 39},{44, 40},{44, 41},{44, 42},{44, 43},
                                     {45, 0},{45, 1},{45, 2},{45, 3},{45, 4},{45, 5},{45, 6},{45, 7},{45, 8},{45, 9},{45, 10},{45, 11},{45, 12},{45, 13},{45, 14},{45, 15},{45, 16},{45, 17},{45, 18},{45, 19},{45, 20},{45, 21},{45, 22},{45, 23},{45, 24},{45, 25},{45, 26},{45, 27},{45, 28},{45, 29},{45, 30},{45, 31},{45, 32},{45, 33},{45, 34},{45, 35},{45, 36},{45, 37},{45, 38},{45, 39},{45, 40},{45, 41},{45, 42},{45, 43},{45, 44},
                                     {46, 0},{46, 1},{46, 2},{46, 3},{46, 4},{46, 5},{46, 6},{46, 7},{46, 8},{46, 9},{46, 10},{46, 11},{46, 12},{46, 13},{46, 14},{46, 15},{46, 16},{46, 17},{46, 18},{46, 19},{46, 20},{46, 21},{46, 22},{46, 23},{46, 24},{46, 25},{46, 26},{46, 27},{46, 28},{46, 29},{46, 30},{46, 31},{46, 32},{46, 33},{46, 34},{46, 35},{46, 36},{46, 37},{46, 38},{46, 39},{46, 40},{46, 41},{46, 42},{46, 43},{46, 44},{46, 45},
                                     {47, 0},{47, 1},{47, 2},{47, 3},{47, 4},{47, 5},{47, 6},{47, 7},{47, 8},{47, 9},{47, 10},{47, 11},{47, 12},{47, 13},{47, 14},{47, 15},{47, 16},{47, 17},{47, 18},{47, 19},{47, 20},{47, 21},{47, 22},{47, 23},{47, 24},{47, 25},{47, 26},{47, 27},{47, 28},{47, 29},{47, 30},{47, 31},{47, 32},{47, 33},{47, 34},{47, 35},{47, 36},{47, 37},{47, 38},{47, 39},{47, 40},{47, 41},{47, 42},{47, 43},{47, 44},{47, 45},{47, 46},
                                     {48, 0},{48, 1},{48, 2},{48, 3},{48, 4},{48, 5},{48, 6},{48, 7},{48, 8},{48, 9},{48, 10},{48, 11},{48, 12},{48, 13},{48, 14},{48, 15},{48, 16},{48, 17},{48, 18},{48, 19},{48, 20},{48, 21},{48, 22},{48, 23},{48, 24},{48, 25},{48, 26},{48, 27},{48, 28},{48, 29},{48, 30},{48, 31},{48, 32},{48, 33},{48, 34},{48, 35},{48, 36},{48, 37},{48, 38},{48, 39},{48, 40},{48, 41},{48, 42},{48, 43},{48, 44},{48, 45},{48, 46},{48, 47},
                                     {49, 0},{49, 1},{49, 2},{49, 3},{49, 4},{49, 5},{49, 6},{49, 7},{49, 8},{49, 9},{49, 10},{49, 11},{49, 12},{49, 13},{49, 14},{49, 15},{49, 16},{49, 17},{49, 18},{49, 19},{49, 20},{49, 21},{49, 22},{49, 23},{49, 24},{49, 25},{49, 26},{49, 27},{49, 28},{49, 29},{49, 30},{49, 31},{49, 32},{49, 33},{49, 34},{49, 35},{49, 36},{49, 37},{49, 38},{49, 39},{49, 40},{49, 41},{49, 42},{49, 43},{49, 44},{49, 45},{49, 46},{49, 47},{49, 48},
                                     {50, 0},{50, 1},{50, 2},{50, 3},{50, 4},{50, 5},{50, 6},{50, 7},{50, 8},{50, 9},{50, 10},{50, 11},{50, 12},{50, 13},{50, 14},{50, 15},{50, 16},{50, 17},{50, 18},{50, 19},{50, 20},{50, 21},{50, 22},{50, 23},{50, 24},{50, 25},{50, 26},{50, 27},{50, 28},{50, 29},{50, 30},{50, 31},{50, 32},{50, 33},{50, 34},{50, 35},{50, 36},{50, 37},{50, 38},{50, 39},{50, 40},{50, 41},{50, 42},{50, 43},{50, 44},{50, 45},{50, 46},{50, 47},{50, 48},{50, 49},
                                     {51, 0},{51, 1},{51, 2},{51, 3},{51, 4},{51, 5},{51, 6},{51, 7},{51, 8},{51, 9},{51, 10},{51, 11},{51, 12},{51, 13},{51, 14},{51, 15},{51, 16},{51, 17},{51, 18},{51, 19},{51, 20},{51, 21},{51, 22},{51, 23},{51, 24},{51, 25},{51, 26},{51, 27},{51, 28},{51, 29},{51, 30},{51, 31},{51, 32},{51, 33},{51, 34},{51, 35},{51, 36},{51, 37},{51, 38},{51, 39},{51, 40},{51, 41},{51, 42},{51, 43},{51, 44},{51, 45},{51, 46},{51, 47},{51, 48},{51, 49},{51, 50},
                                     {52, 0},{52, 1},{52, 2},{52, 3},{52, 4},{52, 5},{52, 6},{52, 7},{52, 8},{52, 9},{52, 10},{52, 11},{52, 12},{52, 13},{52, 14},{52, 15},{52, 16},{52, 17},{52, 18},{52, 19},{52, 20},{52, 21},{52, 22},{52, 23},{52, 24},{52, 25},{52, 26},{52, 27},{52, 28},{52, 29},{52, 30},{52, 31},{52, 32},{52, 33},{52, 34},{52, 35},{52, 36},{52, 37},{52, 38},{52, 39},{52, 40},{52, 41},{52, 42},{52, 43},{52, 44},{52, 45},{52, 46},{52, 47},{52, 48},{52, 49},{52, 50},{52, 51},
                                     {53, 0},{53, 1},{53, 2},{53, 3},{53, 4},{53, 5},{53, 6},{53, 7},{53, 8},{53, 9},{53, 10},{53, 11},{53, 12},{53, 13},{53, 14},{53, 15},{53, 16},{53, 17},{53, 18},{53, 19},{53, 20},{53, 21},{53, 22},{53, 23},{53, 24},{53, 25},{53, 26},{53, 27},{53, 28},{53, 29},{53, 30},{53, 31},{53, 32},{53, 33},{53, 34},{53, 35},{53, 36},{53, 37},{53, 38},{53, 39},{53, 40},{53, 41},{53, 42},{53, 43},{53, 44},{53, 45},{53, 46},{53, 47},{53, 48},{53, 49},{53, 50},{53, 51},{53, 52},
                                     {54, 0},{54, 1},{54, 2},{54, 3},{54, 4},{54, 5},{54, 6},{54, 7},{54, 8},{54, 9},{54, 10},{54, 11},{54, 12},{54, 13},{54, 14},{54, 15},{54, 16},{54, 17},{54, 18},{54, 19},{54, 20},{54, 21},{54, 22},{54, 23},{54, 24},{54, 25},{54, 26},{54, 27},{54, 28},{54, 29},{54, 30},{54, 31},{54, 32},{54, 33},{54, 34},{54, 35},{54, 36},{54, 37},{54, 38},{54, 39},{54, 40},{54, 41},{54, 42},{54, 43},{54, 44},{54, 45},{54, 46},{54, 47},{54, 48},{54, 49},{54, 50},{54, 51},{54, 52},{54, 53},
                                     {55, 0},{55, 1},{55, 2},{55, 3},{55, 4},{55, 5},{55, 6},{55, 7},{55, 8},{55, 9},{55, 10},{55, 11},{55, 12},{55, 13},{55, 14},{55, 15},{55, 16},{55, 17},{55, 18},{55, 19},{55, 20},{55, 21},{55, 22},{55, 23},{55, 24},{55, 25},{55, 26},{55, 27},{55, 28},{55, 29},{55, 30},{55, 31},{55, 32},{55, 33},{55, 34},{55, 35},{55, 36},{55, 37},{55, 38},{55, 39},{55, 40},{55, 41},{55, 42},{55, 43},{55, 44},{55, 45},{55, 46},{55, 47},{55, 48},{55, 49},{55, 50},{55, 51},{55, 52},{55, 53},{55, 54},
                                     {56, 0},{56, 1},{56, 2},{56, 3},{56, 4},{56, 5},{56, 6},{56, 7},{56, 8},{56, 9},{56, 10},{56, 11},{56, 12},{56, 13},{56, 14},{56, 15},{56, 16},{56, 17},{56, 18},{56, 19},{56, 20},{56, 21},{56, 22},{56, 23},{56, 24},{56, 25},{56, 26},{56, 27},{56, 28},{56, 29},{56, 30},{56, 31},{56, 32},{56, 33},{56, 34},{56, 35},{56, 36},{56, 37},{56, 38},{56, 39},{56, 40},{56, 41},{56, 42},{56, 43},{56, 44},{56, 45},{56, 46},{56, 47},{56, 48},{56, 49},{56, 50},{56, 51},{56, 52},{56, 53},{56, 54},{56, 55},
                                     {57, 0},{57, 1},{57, 2},{57, 3},{57, 4},{57, 5},{57, 6},{57, 7},{57, 8},{57, 9},{57, 10},{57, 11},{57, 12},{57, 13},{57, 14},{57, 15},{57, 16},{57, 17},{57, 18},{57, 19},{57, 20},{57, 21},{57, 22},{57, 23},{57, 24},{57, 25},{57, 26},{57, 27},{57, 28},{57, 29},{57, 30},{57, 31},{57, 32},{57, 33},{57, 34},{57, 35},{57, 36},{57, 37},{57, 38},{57, 39},{57, 40},{57, 41},{57, 42},{57, 43},{57, 44},{57, 45},{57, 46},{57, 47},{57, 48},{57, 49},{57, 50},{57, 51},{57, 52},{57, 53},{57, 54},{57, 55},{57, 56},
                                     {58, 0},{58, 1},{58, 2},{58, 3},{58, 4},{58, 5},{58, 6},{58, 7},{58, 8},{58, 9},{58, 10},{58, 11},{58, 12},{58, 13},{58, 14},{58, 15},{58, 16},{58, 17},{58, 18},{58, 19},{58, 20},{58, 21},{58, 22},{58, 23},{58, 24},{58, 25},{58, 26},{58, 27},{58, 28},{58, 29},{58, 30},{58, 31},{58, 32},{58, 33},{58, 34},{58, 35},{58, 36},{58, 37},{58, 38},{58, 39},{58, 40},{58, 41},{58, 42},{58, 43},{58, 44},{58, 45},{58, 46},{58, 47},{58, 48},{58, 49},{58, 50},{58, 51},{58, 52},{58, 53},{58, 54},{58, 55},{58, 56},{58, 57},
                                     {59, 0},{59, 1},{59, 2},{59, 3},{59, 4},{59, 5},{59, 6},{59, 7},{59, 8},{59, 9},{59, 10},{59, 11},{59, 12},{59, 13},{59, 14},{59, 15},{59, 16},{59, 17},{59, 18},{59, 19},{59, 20},{59, 21},{59, 22},{59, 23},{59, 24},{59, 25},{59, 26},{59, 27},{59, 28},{59, 29},{59, 30},{59, 31},{59, 32},{59, 33},{59, 34},{59, 35},{59, 36},{59, 37},{59, 38},{59, 39},{59, 40},{59, 41},{59, 42},{59, 43},{59, 44},{59, 45},{59, 46},{59, 47},{59, 48},{59, 49},{59, 50},{59, 51},{59, 52},{59, 53},{59, 54},{59, 55},{59, 56},{59, 57},{59, 58},
                                     {60, 0},{60, 1},{60, 2},{60, 3},{60, 4},{60, 5},{60, 6},{60, 7},{60, 8},{60, 9},{60, 10},{60, 11},{60, 12},{60, 13},{60, 14},{60, 15},{60, 16},{60, 17},{60, 18},{60, 19},{60, 20},{60, 21},{60, 22},{60, 23},{60, 24},{60, 25},{60, 26},{60, 27},{60, 28},{60, 29},{60, 30},{60, 31},{60, 32},{60, 33},{60, 34},{60, 35},{60, 36},{60, 37},{60, 38},{60, 39},{60, 40},{60, 41},{60, 42},{60, 43},{60, 44},{60, 45},{60, 46},{60, 47},{60, 48},{60, 49},{60, 50},{60, 51},{60, 52},{60, 53},{60, 54},{60, 55},{60, 56},{60, 57},{60, 58},{60, 59},
                                     {61, 0},{61, 1},{61, 2},{61, 3},{61, 4},{61, 5},{61, 6},{61, 7},{61, 8},{61, 9},{61, 10},{61, 11},{61, 12},{61, 13},{61, 14},{61, 15},{61, 16},{61, 17},{61, 18},{61, 19},{61, 20},{61, 21},{61, 22},{61, 23},{61, 24},{61, 25},{61, 26},{61, 27},{61, 28},{61, 29},{61, 30},{61, 31},{61, 32},{61, 33},{61, 34},{61, 35},{61, 36},{61, 37},{61, 38},{61, 39},{61, 40},{61, 41},{61, 42},{61, 43},{61, 44},{61, 45},{61, 46},{61, 47},{61, 48},{61, 49},{61, 50},{61, 51},{61, 52},{61, 53},{61, 54},{61, 55},{61, 56},{61, 57},{61, 58},{61, 59},{61, 60},
                                     {62, 0},{62, 1},{62, 2},{62, 3},{62, 4},{62, 5},{62, 6},{62, 7},{62, 8},{62, 9},{62, 10},{62, 11},{62, 12},{62, 13},{62, 14},{62, 15},{62, 16},{62, 17},{62, 18},{62, 19},{62, 20},{62, 21},{62, 22},{62, 23},{62, 24},{62, 25},{62, 26},{62, 27},{62, 28},{62, 29},{62, 30},{62, 31},{62, 32},{62, 33},{62, 34},{62, 35},{62, 36},{62, 37},{62, 38},{62, 39},{62, 40},{62, 41},{62, 42},{62, 43},{62, 44},{62, 45},{62, 46},{62, 47},{62, 48},{62, 49},{62, 50},{62, 51},{62, 52},{62, 53},{62, 54},{62, 55},{62, 56},{62, 57},{62, 58},{62, 59},{62, 60},{62, 61},
                                     {63, 0},{63, 1},{63, 2},{63, 3},{63, 4},{63, 5},{63, 6},{63, 7},{63, 8},{63, 9},{63, 10},{63, 11},{63, 12},{63, 13},{63, 14},{63, 15},{63, 16},{63, 17},{63, 18},{63, 19},{63, 20},{63, 21},{63, 22},{63, 23},{63, 24},{63, 25},{63, 26},{63, 27},{63, 28},{63, 29},{63, 30},{63, 31},{63, 32},{63, 33},{63, 34},{63, 35},{63, 36},{63, 37},{63, 38},{63, 39},{63, 40},{63, 41},{63, 42},{63, 43},{63, 44},{63, 45},{63, 46},{63, 47},{63, 48},{63, 49},{63, 50},{63, 51},{63, 52},{63, 53},{63, 54},{63, 55},{63, 56},{63, 57},{63, 58},{63, 59},{63, 60},{63, 61},{63, 62}};

//expectation step starts
void normalizing_qmatrix(gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *qmatrix_temp) {
    //set qmatrix_temp & normalize
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
}

void get_eigenvector_and_inverse(gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *qmatrix_temp, gsl_matrix *eigenvector, gsl_matrix *eigenvec_inverse, gsl_vector *eigenvalue) {
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
    gsl_vector_complex_free(eigenvalue_com);

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

void set_matrices(newick_start *start, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, int newick_order_max) {
    std::vector<newick_graph*> next_iterator = start->next;
    gsl_matrix *eigenvector_temp = gsl_matrix_alloc(64, 64);//freed
    gsl_vector *eigenvalue_temp = gsl_vector_alloc(64);//freed
    int newick_order = newick_order_max;
    while (newick_order >= 0) {
        size_t size = next_iterator.size();
        for (size_t num = 0; num < size; num++) {
            if (next_iterator[num]->order == newick_order) {
                gsl_matrix_memcpy(next_iterator[num]->expon_matrix[0],
                                  calculate_expon_matrix(eigenvector, eigen_inverse, eigenvalue, eigenvector_temp, eigenvalue_temp, next_iterator[num]->branch_length));
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

double felsenstein_algorithm(newick_graph *node, char base) {
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
            sum_first += node->previous[0]->felsenstein[num] * gsl_matrix_get(node->previous[0]->expon_matrix[0], base, num);
            sum_second += node->previous[1]->felsenstein[num] * gsl_matrix_get(node->previous[1]->expon_matrix[0], base, num);
        }
        node->felsenstein[base] = sum_first * sum_second;
        return sum_first * sum_second;
    }
}

void conduct_felsenstein(newick_start *start, std::vector<aligned_codon> codon_set, int newick_order_max) {
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
    int newick_order = newick_order_max;
    while (newick_order >= 0) {
        int size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            if (next_iteration[num]->order == newick_order) {
                for (char codon = 0; codon < 64; codon++) {
                    felsenstein_algorithm(next_iteration[num], codon);
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
                        fel_upper += next_iteration[num]->next->upper[base_first]
                                     * gsl_matrix_get(next_iteration[num]->expon_matrix[0], base_first, base)
                                     * next_iteration[num]->next->previous[(next_iteration[num] == next_iteration[num]->next->previous[0]) ? 1 : 0]->felsenstein[base_second]
                                     * gsl_matrix_get(next_iteration[num]->next->previous[(next_iteration[num] == next_iteration[num]->next->previous[0]) ? 1 : 0]->expon_matrix[0], base_first, base_second);
                    }
                }
                next_iteration[num]->upper[base] = fel_upper;
            }
            if (!(next_iteration[num]->previous.empty())) {
                next_iteration.push_back(next_iteration[num]->previous[0]);
                next_iteration.push_back(next_iteration[num]->previous[1]);
            }
        }
        next_iteration.erase(next_iteration.begin(), next_iteration.begin() + size);
    }
}

void update_upper(newick_start *start, double *codon_freq, int newick_order_max) {
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

    int newick_order = newick_order_max;
    //update upper
    while (newick_order >= 0) {
        int size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            if (next_iteration[num]->order == newick_order) {
                for (char base = 0; base < 64; base++) {
                    next_iteration[num]->upper[base] = next_iteration[num]->felsenstein[base] * next_iteration[num]->upper[base] / denominator;
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

void calculate_expectation(newick_start *start, int newick_order_max) {
    std::vector<newick_graph*> next_iteration = start->next;
    int newick_order = newick_order_max;
    while (newick_order >= 0) {
        int size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            if (next_iteration[num]->order == newick_order) {
                for (char base_next = 0; base_next < 64; base_next++) {
                    double denominator = 0.0;
                    for (char base_acc = 0; base_acc < 64; base_acc++) {
                        denominator += next_iteration[num]->felsenstein[base_acc] *
                                       gsl_matrix_get(next_iteration[num]->expon_matrix[0], base_next, base_acc);
                    }
                    //Todo: if denominator = 0, for sure numerator is 0, so make the expectation 0 (no need to calculate) (check if this is the right solution later)
                    if (denominator == 0) {
                        for (char base_curr = 0; base_curr < 64; base_curr++) {
                            next_iteration[num]->expectation[64 * base_next + base_curr] = 0.0;
                        }
                        continue;
                    }
                    for (char base_curr = 0; base_curr < 64; base_curr++) {
                        next_iteration[num]->expectation[64 * base_next + base_curr] +=
                                next_iteration[num]->felsenstein[base_curr] *
                                gsl_matrix_get(next_iteration[num]->expon_matrix[0], base_next, base_curr) *
                                next_iteration[num]->next->upper[base_next] / denominator;
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

//Integrating all the expectation steps
void conduct_expectation_step(std::vector<std::vector<aligned_codon>> aligned_codon_set, newick_start *start, newick_graph *end,
                              gsl_matrix *qmatrix, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, double *codon_freq, int newick_order_max) {
    std::vector<newick_graph*> next_iterator = start->next;
    for (size_t num = 0; num < next_iterator.size(); num++) {
        for (char base = 0; base < 64; base++) {
            next_iterator[num]->base[base] = true;
        }
    }

    int newick_order = newick_order_max;
    while (newick_order >= 0) {
        size_t size = next_iterator.size();
        for (size_t num = 0; num < size; num++) {
            if (next_iterator[num]->order == newick_order) {
                for (char base_first = 0; base_first < 64; base_first++) {
                    for (char base_second = 0; base_second < 64; base_second++) {
                        next_iterator[num]->expectation[64 * base_first + base_second] = 0.0;
                    }
                }
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

    gsl_matrix *qmatrix_temp = gsl_matrix_alloc(64, 64);//freed

    normalizing_qmatrix(qmatrix, codon_freq, qmatrix_temp);

    get_eigenvector_and_inverse(qmatrix, codon_freq, qmatrix_temp, eigenvector, eigen_inverse, eigenvalue);

    set_matrices(start, eigenvector, eigen_inverse, eigenvalue, newick_order_max);

    for (size_t num = 0; num < aligned_codon_set.size(); num++) {
        std::cout << num << std::endl;
        conduct_felsenstein(start, aligned_codon_set[num], newick_order_max);

        calculate_upper(end, codon_freq);

        update_upper(start, codon_freq, newick_order_max);

        calculate_expectation(start, newick_order_max);
    }
}

//copying maximization_step.cpp
void calculate_derivative(gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, newick_start *start, gsl_matrix *gradient, int newick_order_max) {
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
        thread_local int newick_order = newick_order_max;
        thread_local std::vector<newick_graph*> next_iteration = start->next;
        thread_local gsl_matrix *F = gsl_matrix_alloc(64, 64);//freed
        thread_local size_t row_tar = num_to_coordinate[diff][0];
        thread_local size_t col_tar = num_to_coordinate[diff][1];
        while (newick_order >= 0) {
            size_t size = next_iteration.size();
            for (size_t num = 0; num < size; num++) {
                if (next_iteration[num]->order == newick_order) {
                    gsl_matrix_set_all(dxepon_matrix[diff], 0.0);
                    for (size_t row = 0; row < 64; row++) {
                        for (size_t col = 0; col < 64; col++) {
                            if (row == col) {
                                gsl_matrix_set(dxepon_matrix[diff], row, col,
                                               (2 * codon_freq[row_tar] * codon_freq[col_tar] * diagonal[row] -
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
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigen_inverse, dxepon_matrix[diff], 0.0, dxepon_matrix[diff]);
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dxepon_matrix[diff], eigenvector, 0.0, dxepon_matrix[diff]);
                    for (size_t row = 0; row < 64; row++) {
                        for (size_t col = 0; col < 64; col++) {
                            if (abs(gsl_vector_get(eigenvalue, row) - gsl_vector_get(eigenvalue, col)) < 1.0e-14) {
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
                    gsl_matrix_free(F);
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigenvector, dxepon_matrix[diff], 0.0, dxepon_matrix[diff]);
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dxepon_matrix[diff], eigen_inverse, 0.0, dxepon_matrix[num]);
                    //putting negative of gradient
                    for (size_t row = 0; row < 64; row++) {
                        for (size_t col = 0; col < 64; col++) {
                            gsl_matrix_set(gradient, diff, 0, gsl_matrix_get(gradient, diff, 0) - next_iteration[num]->expectation[64 * row + col]
                                                              * gsl_matrix_get(dxepon_matrix[diff], row, col) / gsl_matrix_get(next_iteration[num]->expon_matrix[0], row, col));
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

double calculate_maximization_function(newick_start *start, short int index, int newick_order_max) {
    double result = 0;
    std::vector<newick_graph*> next_iteration = start->next;
    int newick_order = newick_order_max;
    while(newick_order >= 0) {
        size_t size = next_iteration.size();
        for (size_t num = 0; num < size; num++) {
            if (next_iteration[num]->order == newick_order) {
                for (size_t row = 0; row < 64; row++) {
                    for (size_t col = 0; col < 64; col++) {
                        result -= next_iteration[num]->expectation[64 * row + col] *
                                  log(gsl_matrix_get(next_iteration[num]->expon_matrix[index], row, col));
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
    return result;
}

void update_expon_matrix(newick_start *start, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, int newick_order_max) {
    std::vector<newick_graph*> next_iterator = start->next;
    gsl_matrix *eigenvector_temp = gsl_matrix_alloc(64, 64);//freed
    gsl_vector *eigenvalue_temp = gsl_vector_alloc(64);//freed
    int newick_order = newick_order_max;
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

//Todo: do I have to limit the length of each step? If I do, then how much?
double backtracking(newick_start *start, gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *direction, gsl_matrix *gradient, double &lamb, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, int newick_order_max) {
    gsl_matrix *qmatrix_temp = gsl_matrix_alloc(64, 64);//freed
    gsl_matrix_memcpy(qmatrix_temp, qmatrix);
    //TODO: Do I have to normalize in this step??-->seems like has to be normalized (equation needs normalization)
    double funcg_ori = calculate_maximization_function(start, 0, newick_order_max);
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
        row = num_to_coordinate[num][0];
        col = num_to_coordinate[num][1];
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
        update_expon_matrix(start, eigenvector, eigen_inverse, eigenvalue, newick_order_max);

        //calculate g(lambda)
        funcg = calculate_maximization_function(start, 1, newick_order_max);
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

double quasi_Newton_method(newick_start *start, gsl_matrix *qmatrix, double *codon_freq, gsl_matrix *eigenvector, gsl_matrix *eigen_inverse, gsl_vector *eigenvalue, int newick_order_max) {
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
    calculate_derivative(qmatrix, codon_freq, eigenvector, eigen_inverse, eigenvalue, start, gradient, newick_order_max);
    gsl_matrix *dx = gsl_matrix_alloc(2016, 1);//freed
    gsl_matrix *x = gsl_matrix_alloc(2016, 1);//freed
    gsl_matrix *x_new = gsl_matrix_alloc(2016, 1);//freed
    for (int num = 0; num < 2016; num++) {
        gsl_matrix_set(dx, num, 0, -gsl_matrix_get(gradient, num, 0));
    }
    int row_map, col_map;
    for (int num = 0; num < 2016; num++) {
        row_map = num_to_coordinate[num][0];
        col_map = num_to_coordinate[num][1];
        gsl_matrix_set(x, num, 0, gsl_matrix_get(qmatrix, row_map, col_map));
    }

    gsl_matrix *qmatrix_new = gsl_matrix_alloc(64, 64);//freed
    double lambda = 1.0;
    double function_value;
    //start loop
    for (int its = 0; its < iterate_max; its++) {
        function_value = backtracking(start, qmatrix, codon_freq, dx, gradient, lambda, eigenvector, eigen_inverse, eigenvalue, newick_order_max);
        for (int num = 0; num < 2016; num++) {
            gsl_matrix_set(dx, num, 0, lambda * gsl_matrix_get(dx, num, 0));
            row_map = num_to_coordinate[num][0];
            col_map = num_to_coordinate[num][1];
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
            //change all the negative entries into non-negative entries (Israel, Rosenthal & Wei (2001))
            double g_i = 0.0;
            double b_i = 0.0;
            for (size_t row = 0; row < 64; row++) {
                for (size_t col = 0; col < 64; col++) {
                    if (row == col) {
                        g_i += abs(gsl_matrix_get(qmatrix, row, col));
                    } else if (gsl_matrix_get(qmatrix, row, col) >= 0) {
                        g_i += gsl_matrix_get(qmatrix, row, col);
                    } else {
                        b_i -= gsl_matrix_get(qmatrix, row, col);
                    }
                }
                for (size_t col = 0; col < 64; col++) {
                    if (row != col && gsl_matrix_get(qmatrix, row, col) < 0) {
                        gsl_matrix_set(qmatrix, row, col, 0.0);
                    } else if (g_i > 0) {
                        gsl_matrix_set(qmatrix, row, col, gsl_matrix_get(qmatrix, row, col) - b_i * abs(gsl_matrix_get(qmatrix, row, col)) / g_i);
                    }
                }
                g_i = 0.0;
                b_i = 0.0;
            }
            return function_value;
        }
        calculate_derivative(qmatrix_new, codon_freq, eigenvector, eigen_inverse, eigenvalue, start, gradient_new, newick_order_max);
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