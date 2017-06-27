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
#include <string>
#include "mpi.h"
//#include "nodeTree.h"
#include "controlStruct.h"
#include "csrSpMV.h"

void csrClusterSplit_Overflow(char *filename, bool colMajor, std::string distributionMethod,
                              std::vector<int>& origin_row, std::vector<int>& origin_col,
                              std::vector<double>& origin_data, std::vector<std::vector<int> >& temp_row,
                              std::vector<std::vector<int> >& temp_col, std::vector<std::vector<double> >& temp_data,
                              std::vector<int>& colMasterTemp_row, std::vector<int>& colMasterTemp_col,
                              std::vector<double>& colMasterTemp_data, int& rowCount, int& colCount, int& nonZeros,
                              int& colsPerNode, int clusterRows, int clusterCols);

void csrClusterSplit_SplitMatrix(controlData controlData, std::vector<csrSpMV*> clusterColData);

void csrClusterSplit_ElementBalanced(char *filename, bool colMajor, std::string distributionMethod, int processCount,
                                     std::vector<int>& origin_row, std::vector<int>& origin_col,
                                     std::vector<double>& origin_data, std::vector<std::vector<int> >& temp_row,
                                     std::vector<std::vector<int> >& temp_col,
                                     std::vector<std::vector<double> >& temp_data,
                                     int& rowCount, int& colCount, int& nonZeros, int& colsPerNode, int clusterRows,
                                     int clusterCols, std::vector<std::vector <std::vector <int> > >& nodeRowOwnership);

void balanceDistribution(int processCount, int nodeBalanceElementCount, std::vector<std::vector <int> >& rowLengths,
                         std::vector<std::vector <std::vector <int> > >& nodeRowOwnership);

void splitDenseVector_ElementBalanced(std::vector<double> denseVector,
                                      std::vector<std::vector<double> >& splitDenseVector,
                                      std::vector<std::vector <std::vector <int> > >& nodeRowOwnership);

#endif //DISTRUBUTED_SPMV_CSRCLUSTERSPLIT_H
