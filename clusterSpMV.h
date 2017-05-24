//
// Created by brianpage on 5/15/17.
//
#include <vector>
#include <omp.h>
#include <iostream>

#ifndef DISTRUBUTED_SPMV_CLUSTERSPMV_H
#define DISTRUBUTED_SPMV_CLUSTERSPMV_H

void clusterSpMV(int ompThreads, std::vector<int> csr_row, std::vector<int> csr_col, std::vector<double> csr_data,
                 double *denseVector, double *nodeResult, int rowsPerNode, bool colMajor);

#endif //DISTRUBUTED_SPMV_CLUSTERSPMV_H
