#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cmath>
#include "mpi.h"
#include "mm_to_csr.h"

int main(int argc, char *argv[]) {

    //Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int processCount;
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    // Get the rank of the process
    int myId;
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Verify MPI initialised properly
    if (processCount == 1) {
        printf("MPI Initialization failed -- KILLING!");
        MPI_Finalize();
        exit(0);
    }

    //**************************************************//
    //  Create comm for each column of compute nodes    //
    //**************************************************//
    // get the number of MPI processes for work split and distribution
    int clusterRows = sqrt(processCount);
    int myCol = myId % clusterRows;

    // Create a duplicate of MPI_COMM_WORLD that can be used to split
    MPI_Comm dupCommWorldCol;
    MPI_Comm_dup(MPI_COMM_WORLD, &dupCommWorldCol);

    // Split dupCommWorld into comm's for each node column, based on each processes original myId
    MPI_Comm col_comm;
    MPI_Comm_split(dupCommWorldCol, myCol, myId, &col_comm);

    int col_rank, col_size;
    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_size(col_comm, &col_size);

    printf("WORLD RANK/SIZE: %d/%d \t COL RANK/SIZE: %d/%d\n", myId, processCount, col_rank, col_size);
    //*************************************//
    //     End creation of column comms    //
    //************************************ //

    //*********************************************//
    //  Create comm for each row of compute nodes  //
    //*********************************************//
    // get the number of MPI processes for work split and distribution
    int myRow = myId / clusterRows;
    // Create a duplicate of MPI_COMM_WORLD that can be used to split
    MPI_Comm dupCommWorldRow;
    MPI_Comm_dup(MPI_COMM_WORLD, &dupCommWorldRow);

    // Split dupCommWorld into comm's for each node column, based on each processes original myId
    MPI_Comm row_comm;
    MPI_Comm_split(dupCommWorldRow, myRow, myId, &row_comm);

    int row_rank, row_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);

    printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n", myId, processCount, row_rank, row_size);
    //**********************************//
    //     End creation of row comms    //
    //**********************************//

    // create vectors to hold sparse matrix data once converted to CSR format
    std::vector<int> csr_row;
    std::vector<int> csr_col;
    std::vector<double> csr_data;
    int rowCount, colCount, nonZeros;

    if (!myId) {
        //convert sparse matrix from Matrix Market to Compressed Sparse Row format
        MMCOO_to_CSR(argv[1], csr_row, csr_col, csr_data, rowCount, colCount, nonZeros);
    }

    // Let everyone know the size of the vectors we will be sending them
    MPI_Bcast(&rowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&colCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nonZeros, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // resize vectors on non master nodes, as they have not yet done this.
    if (myId) {
        csr_row.resize(rowCount);
        csr_col.resize(colCount);
        csr_data.resize(nonZeros);
    }

    //if (!myId) std::cout << "Broadcasting array data from process myId to all other processes" << std::endl;
    //send data to every cluster node
    MPI_Bcast(&csr_row[0], rowCount, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&csr_col[0], colCount, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&csr_data[0], nonZeros, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    int vecElementsPerCol = rowCount / clusterRows;
    int matrixRowsPerNode = rowCount / clusterRows;
    double *denseVector;
    if (!myId) {
        denseVector = new double[rowCount];

        for (int i = 0; i < rowCount; i++) {
            denseVector[i] = rand();
        }
    }

    double denseVecColPiece[vecElementsPerCol];
    //std::cout << "denseVector size = " << rowCount << ", denseVecColPeice size = " << rowsPerCol << std::endl;
    if (myId < clusterRows) {
        //std::cout << "sending col vector portion to column roots" << std::endl;
        MPI_Scatter(denseVector, vecElementsPerCol, MPI_DOUBLE, &denseVecColPiece, vecElementsPerCol, MPI_DOUBLE, 0, row_comm);
        MPI_Barrier(row_comm);
    }

    //if (!myId) std::cout << "broadcasting col vector portion to all column nodes" << std::endl;
    MPI_Bcast(&denseVecColPiece, vecElementsPerCol, MPI_DOUBLE, 0, col_comm);
    MPI_Barrier(col_comm);
    //if (!myId) std::cout << "done sending all data " << std::endl;

    // Do Sparse Matrix Dense Vector Multiplication on all nodes
    int myWorkRowStart = myRow * matrixRowsPerNode;
    int myWorkColStart = myCol * vecElementsPerCol;

    if (!myId)std::cout << "myId = " << myId << ", rows = "<<  myWorkRowStart << "-" << myWorkRowStart+matrixRowsPerNode << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    if (!myId) std::cout << "myId = " << myId << ", col = "<<  myWorkColStart << "-" << myWorkColStart+vecElementsPerCol << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);


    int *myCSR_Row;
    int *myCSR_Col;
    double *myCSR_Data;
/*
    if (!myId) std::cout << csr_row[0] << std::endl;
    if (!myId) std::cout << csr_col[0] << std::endl;
    if (!myId) std::cout << csr_data[0] << std::endl << std::endl << std::endl;

    if (!myId) std::cout << csr_row[1] << std::endl;
    if (!myId) std::cout << csr_col[1] << std::endl;
    if (!myId) std::cout << csr_data[1] << std::endl;

    for (int i = 0; i < 2; i++){
        if (!myId) std::cout << csr_row[i] << std::endl;
    }
*/

    myCSR_Row = new int[matrixRowsPerNode];
    int firstData = csr_row[myWorkRowStart];
    int lastData = csr_row[myWorkRowStart+matrixRowsPerNode]-1;
    int numOfElements = lastData-firstData;
    std::cout << firstData << "-" << lastData << ", numOfElements = " << numOfElements << std::endl;
    myCSR_Col = new int[numOfElements];
    myCSR_Data = new double[numOfElements];

    int count = 0;
    for (int i = myWorkRowStart; i < myWorkRowStart+matrixRowsPerNode; i++){
        myCSR_Row[count] = csr_row[i];
        count++;
    }

    count = 0;
    for (int i = firstData; i < numOfElements; i++){
        //std::cout << count << ", " << i << std::endl;
        myCSR_Col[count] = csr_col[i];
        myCSR_Data[count] = csr_data[i];
        count++;
    }

    int k = 0;
    double colY[matrixRowsPerNode] = {0.0};
    for (int i = 0; i < matrixRowsPerNode; i++){
        int j = 0;

        while ((myCSR_Col[myCSR_Row[i]+j] >= myWorkColStart) &&
                (myCSR_Row[myCSR_Row[i]+j] <= myWorkColStart+vecElementsPerCol)){
            std::cout << "colY[" << k << "] += myCSR_Data[" << myCSR_Row[i] << " + "
                      << j << "] * denseVecColPiece[" << j << "]" << std::endl;
            colY[k] += myCSR_Data[myCSR_Row[i] + j] * denseVecColPiece[j];
            j++;
        }
        k++;
    }


    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "DONE CALCULATING SHIT" << std::endl;
/*
    double *rowY;
    double *rowYSum;
    if ((myId % clusterRows) == 0){
        rowY = new double[matrixRowsPerNode*clusterRows];
        MPI_Gather(&rowY, matrixRowsPerNode, MPI_DOUBLE, colY, matrixRowsPerNode, MPI_DOUBLE, 0, row_comm);

        rowYSum = new double[matrixRowsPerNode];
        double temp = 0.0;
        for (int i = 0; i < matrixRowsPerNode; i++){
            for (int j = 0; j < matrixRowsPerNode*clusterRows; j = j+matrixRowsPerNode){
                temp += rowY[j];
            }
            rowYSum[i] = temp;
        }
    }
*/

    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&row_comm);
    MPI_Finalize();

    return 0;
}