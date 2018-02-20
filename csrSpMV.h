//
// Created by brianpage on 6/22/17.
//

#ifndef DISTRUBUTED_SPMV_CSRSPMV_H
#define DISTRUBUTED_SPMV_CSRSPMV_H

#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include "controlStruct.h"
#include <chrono>
#include <thread>
#include <omp.h>

class csrSpMV {
    public:

		std::vector <int> processData; // 0 is nnz per row, 1 is # of rows for given process, 2 is densevec length
		std::vector <int> assignedRowIds;
		std::vector <std::vector <int> > assignedRowIds2d;
		std::vector <int> csrRows;
		std::vector <int> csrCols;
		std::vector <double> csrData;
		std::vector <double> denseVec;
        std::vector <double> result;

        void nodeSpMV(controlData control, std::vector <double>& result);
        void masterOnlySpMV(controlData control, std::vector<int>& seqDist);
		csrSpMV();                          // generic constructor
		csrSpMV(const csrSpMV& objToCopy);  //copy constructor
        ~csrSpMV();                         // destructor
		void rebase(int colAdjustment);
		void rebase_balanced();
};

#endif //DISTRUBUTED_SPMV_CSRSPMV_H
