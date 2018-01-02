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
		std::vector <int> csrRows;
		std::vector <int> csrCols;
		std::vector <double> csrData;
		std::vector <double> denseVec;
        std::vector <double> result;

        void nodeSpMV(controlData control, std::vector <double>& result);
        void masterOnlySpMV(controlData control);
		csrSpMV();                          // generic constructor
		csrSpMV(const csrSpMV& objToCopy);  //copy constructor
        ~csrSpMV();                         // destructor
		void rebase(int colAdjustment);
};

#endif //DISTRUBUTED_SPMV_CSRSPMV_H
