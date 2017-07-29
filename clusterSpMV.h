//
// Created by brianpage on 5/15/17.
//
#include <vector>
#include <omp.h>
#include <iostream>

#ifndef DISTRUBUTED_SPMV_CLUSTERSPMV_H
#define DISTRUBUTED_SPMV_CLUSTERSPMV_H

// Not currently implemented in this branch
void clusterSpMV(int ompThreads, std::vector<int> csr_row, std::vector<int> csr_col, std::vector<double> csr_data,
                 double *denseVector, double *nodeResult, int rowsPerNode, bool colMajor);

void clusterSpMV_SplitMatrix(int ompThreads, std::vector<int> csr_row, std::vector<int> csr_col,
                             std::vector<double> csr_data, std::vector<double> denseVector,
                             std::vector<double>& nodeResult, bool colMajor);

// Not currently implemented in this branch
void clusterSpMV_ElementBalanced(int ompThreads, std::vector<int> csr_row, std::vector<int> csr_col,
                                 std::vector<double> csr_data, std::vector<double> denseVector,
                                 std::vector<double>& nodeResult, bool colMajor);

#endif //DISTRUBUTED_SPMV_CLUSTERSPMV_H
