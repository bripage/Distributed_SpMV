//
// Created by brianpage on 5/15/17.
//
#ifndef DISTRUBUTION_H
#define DISTRUBUTION_H

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
#include "controlStruct.h"
#include "csrSpMV.h"
#include <ctime>
#include <cstdlib>

struct Element {
	int row;
	int col;
	double data;

	Element(int r, int c, double d) : row(r), col(c), data(d) {}
};

struct row {
	int rowLength, rowId, processAssignment;
	std::vector <int> rowIds;
	std::vector <int> colIds;
	std::vector <double> data;
};

void distribution_SplitMatrix(controlData& controlData, std::vector<csrSpMV*>& clusterColData);
void distribution_Balanced(controlData& control, std::vector<csrSpMV*>& clusterColData);

#endif //DISTRUBUTION_H
