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
    std::vector<int> csr_row, masterRow;
    std::vector<int> csr_col, masterCol;
    std::vector<double> csr_data, masterData;
    int rowCount, colCount, nonZeros;

    if (!myId) {
        //convert sparse matrix from Matrix Market to Compressed Sparse Row format
        MMCOO_to_CSR(argv[1], csr_row, csr_col, csr_data, rowCount, colCount, nonZeros);
    }

    // Let everyone know the size of the vectors we will be sending them
    //MPI_Bcast(&rowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Bcast(&colCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Bcast(&nonZeros, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // create arrays for individual nodes to hold their data.
    int nodeColCount;
    int rowsPerNode = rowCount / clusterRows;
    if (!myId) {
        std::vector<int> colToSend;
        std::vector<double> dataToSend;
        int rowToSend[rowsPerNode];


        for (int i = 1; i < processCount; i++) {
            int count = 0;
            for (int j = i*rowsPerNode; j < (i+1)*rowsPerNode; j++){
                for (int k = 0; k < rowsPerNode; k++){
                    colToSend.push_back(csr_col[csr_row[j]+k]);
                    dataToSend.push_back(csr_data[csr_row[j]+k]);
                }
                rowToSend[count] = csr_row[j] - csr_row[i*rowsPerNode];
                count++;
            }

            nodeColCount = colToSend.size();

            MPI_Send(&nodeColCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&rowToSend[0], rowsPerNode, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&colToSend[0], nodeColCount, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&dataToSend[0], nodeColCount, MPI_INT, i, 0, MPI_COMM_WORLD);

            colToSend.clear();
            dataToSend.clear();
        }

        int count = 0;
        for (int j = 0; j < rowsPerNode; j++){
            for (int k = 0; k < rowsPerNode; k++){
                masterCol.push_back(csr_col[csr_row[j]+k]);
                masterData.push_back(csr_data[csr_row[j]+k]);
            }
            masterRow.push_back(csr_row[j] - csr_row[0]);
            count++;
        }

        csr_col.clear();
        csr_row.clear();
        csr_data.clear();
    }

    if(myId) {
        MPI_Recv(&nodeColCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        csr_col.resize(nodeColCount);
        csr_data.resize(nodeColCount);

        MPI_Recv(&csr_col, nodeColCount, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&csr_data, nodeColCount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }




    //if (!myId) std::cout << "Broadcasting array data from process myId to all other processes" << std::endl;
    //send data to every cluster node
    /*
    MPI_Bcast(&csr_row[0], rowCount, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&csr_col[0], colCount, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&csr_data[0], nonZeros, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    */
    int vecElementsPerCol = rowCount / clusterRows;
    int matrixRowsPerNode = rowCount / clusterRows;
    double denseVector[rowCount];
    if (!myId) {
        //denseVector = new double[rowCount];

        for (int i = 0; i < rowCount; i++) {
            denseVector[i] = rand();
        }
    }

    /*
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

    std::cout << "myId = " << myId << ", rows = "<<  myWorkRowStart << "-" << myWorkRowStart+matrixRowsPerNode << std::endl;
    std::cout << "myId = " << myId << ", col = "<<  myWorkColStart << "-" << myWorkColStart+vecElementsPerCol << std::endl;
*/
    MPI_Bcast(&denseVector[0], rowCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    double nodeResult[rowsPerNode] = {0.0};

    if(!myId){
        for (int i = 0; i < rowsPerNode; i++){
            for (int j = 0; j < nodeColCount; j++){
                nodeResult[i] += masterData[j] * denseVector[masterCol[j]];
            }
        }
    } else {
        for (int i = 0; i < rowsPerNode; i++){
            for (int j = 0; j < nodeColCount; j++){
                nodeResult[i] += csr_data[j] * denseVector[csr_col[j]];
            }
        }
    }

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