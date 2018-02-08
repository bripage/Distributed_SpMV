//
// Created by brianpage on 6/22/17.
//

#ifndef DISTRUBUTED_SPMV_CONTROLSTRUCT_H
#define DISTRUBUTED_SPMV_CONTROLSTRUCT_H

#include <string>
#include <mpi.h>

struct controlData {
    char *matrixFile;
    char *vectorFile;
    int rowCount;   // from the MMF file
    int colCount;   // from the MMF file
    int nonZeros;   // from the MMF file
    bool masterOnly = false;
    bool colMajor = false;  // column or row major order?
	bool barrier = false;
    bool verify = false;
	bool debug = false;
	bool matrixInfo = false;
	bool distributionInfo = false;
    int ompThreads = 1; // number of OpenMP threads to run on each MPI process (cluster node)
    int distributionMethod; // colBalanced, splitMatrix, overflow
    int rowsPerNode;
    int colsPerNode;
    int clusterRows;
    int clusterCols;
    int lastClusterRowRowStart; // can be calculated locally
    int lastClusterColColStart; // can be calculated locally
    int elementCount;
	int maxRowsAssigned;

    //MPI Stuff
    int processCount;
    int myId;
    int row_rank, row_size, myRow;
    int col_rank, col_size, myCol;
    MPI_Comm row_comm;
    MPI_Comm col_comm;
};

#endif //DISTRUBUTED_SPMV_CONTROLSTRUCT_H
