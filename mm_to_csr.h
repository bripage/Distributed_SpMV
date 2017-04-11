//
// Created by brianpage on 4/6/17.
//
#ifndef MM_TO_CSR
#define MM_TO_CSR

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>

void MMCOO_to_CSR(char *filename, std::vector<int> &csr_row, std::vector<int> &csr_col, std::vector<double> &csr_data,
                  int &rowCount, int &colCount, int &nonZeros);
#endif MM_TO_CSR