//
// Created by brianpage on 5/15/17.
//
#ifndef DISTRUBUTED_SPMV_CSRCLUSTERSPLIT_H
#define DISTRUBUTED_SPMV_CSRCLUSTERSPLIT_H

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <omp.h>
#include "mpi.h"
#include "mm_to_csr.h"
#include <string>


void csrClusterSplit(char *filename, bool colMajor, std::string distributionMethod, std::vector<int>& origin_row,
                     std::vector<int>& origin_col, std::vector<double>& origin_data,
                     std::vector<std::vector<int> >& temp_row,
                     std::vector<std::vector<int> >& temp_col, std::vector<std::vector<double> >& temp_data,
                     std::vector<int>& colMasterTemp_row, std::vector<int>& colMasterTemp_col,
                     std::vector<double>& colMasterTemp_data, int& rowCount, int& colCount, int& nonZeros,
                     int& colsPerNode, int clusterRows, int clusterCols);

void overflowSplit(std::vector<int>& origin_row, std::vector<int>& origin_col,  std::vector<double>& origin_data,
                   std::vector<std::vector<int> >& temp_row, std::vector<std::vector<int> >& temp_col,
                   std::vector<std::vector<double> >& temp_data, std::vector<int>& colMasterTemp_row,
                   std::vector<int>& colMasterTemp_col, std::vector<double>& colMasterTemp_data, int& rowCount,
                   int& colsPerNode, int clusterRows, int clusterCols);

void columnBalancedSplit(std::vector<int>& origin_row, std::vector<int>& origin_col,  std::vector<double>& origin_data,
                   std::vector<std::vector<int> >& temp_row, std::vector<std::vector<int> >& temp_col,
                   std::vector<std::vector<double> >& temp_data, std::vector<int>& colMasterTemp_row,
                   std::vector<int>& colMasterTemp_col, std::vector<double>& colMasterTemp_data, int& rowCount,
                   int& colsPerNode, int clusterRows, int clusterCols);

#endif //DISTRUBUTED_SPMV_CSRCLUSTERSPLIT_H
