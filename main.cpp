#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include "csrClusterSplit.h"
#include "clusterSpMV.h"

int main(int argc, char *argv[]) {

    int ompThreads = atoi(argv[2]);
    std::string argTemp = argv[3];
    bool masterOnly;
    if (argTemp == "true"){
        masterOnly = true;
    } else {
        masterOnly = false;
    }

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

    if (processCount == 1 && masterOnly == false) {
        printf("MPI Initialization failed -- KILLING!");
        MPI_Finalize();
        exit(0);
    }

        //**************************************************//
        //  Create comm for each column of compute nodes    //
        //**************************************************//
        // get the number of MPI processes for work split and distribution
        int clusterRows = sqrt(processCount);
        int clusterCols = clusterRows; // works becuase cluster is expected to be square
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

        //printf("WORLD RANK/SIZE: %d/%d \t COL RANK/SIZE: %d/%d\n", myId, processCount, col_rank, col_size);
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

        //printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n", myId, processCount, row_rank, row_size);
        //**********************************//
        //     End creation of row comms    //
        //**********************************//

    // create vectors to hold sparse matrix data once converted to CSR format
    std::vector<int> origin_row, csr_row, colMasterTemp_row, origin_col, csr_col, colMasterTemp_col;
    std::vector<double> origin_data, csr_data, colMasterTemp_data;
    std::vector<std::vector<int> > temp_row,temp_col;
    std::vector<std::vector<double> > temp_data;


    int rowCount, colCount, nonZeros;
    int colsPerNode;

    if(!myId) {
        csrClusterSplit(argv[1], origin_row, origin_col, origin_data, temp_row, temp_col, temp_data, colMasterTemp_row,
                        colMasterTemp_col, colMasterTemp_data, rowCount, colCount, nonZeros, colsPerNode, clusterRows,
                        clusterCols);

        //std::cout << "rowCount = " << rowCount << ", nonZeros = " << nonZeros  << std::endl;
        //std::cout << "origin_row.size() = " << origin_row.size() << ", origin_col.size() = " << origin_col.size() << ", origin_data.size() = " << origin_data.size() << std::endl;
        //std::cout << "temp_row.size() = " << temp_row.size() << ", temp_col.size() = " << temp_col.size() << ", temp_data.size() = " << temp_data.size() << std::endl;

    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&rowCount, 1, MPI_INT, 0, MPI_COMM_WORLD); // need this in order to allocate denseVector memory
    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << myId << "creating Dense Vector" << std::endl;
    double denseVector[rowCount];
    if (!myId) {
        std::cout << "Populating Dense Vector" << std::endl;
        for (int i = 0; i < rowCount; i++) {
            denseVector[i] = 1.0;
        }
    }


    double *SpMVResult, *masterOnly_SpMV, temp;
    int elementsPerRow = rowCount/clusterRows;
    int numToGather = elementsPerRow * clusterRows;

    // perform SpMV on master only, without any segmentation or distribution
    // for verification purposes.
    if(!myId) {
        std::cout << "Performing Master Only SpMV" << std::endl;
        masterOnly_SpMV = new double[rowCount];
        for (int i = 0; i < rowCount; i++) {
            temp = 0.0;
            if (i != rowCount - 1) {
                for (int j = origin_row[i]; j < origin_row[i + 1]; j++) {
                    temp += origin_data[j] * denseVector[origin_col[j]];
                }
            } else {
                for (int j = origin_row[i]; j < origin_data.size(); j++) {
                    temp += origin_data[j] * denseVector[origin_col[j]];
                }
            }
            masterOnly_SpMV[i] = temp;
        }

        origin_row.clear();
        origin_col.clear();
        origin_data.clear();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // create arrays for individual nodes to hold their data.
    int nodeColCount, colRowCount;
    int rowsPerNode = rowCount / clusterRows;

    if(!myId) {
        // Send column data to column masters
        for (int i = 1; i < clusterRows; i++) {
            nodeColCount = temp_data[i].size();
            colRowCount = temp_row[i].size();

            MPI_Send(&nodeColCount, 1, MPI_INT, i, 0, row_comm);
            MPI_Send(&colRowCount, 1, MPI_INT, i, 0, row_comm);
            MPI_Send(&temp_row[i][0], rowCount, MPI_INT, i, 0, row_comm);
            MPI_Send(&temp_col[i][0], nodeColCount, MPI_INT, i, 0, row_comm);
            MPI_Send(&temp_data[i][0], nodeColCount, MPI_DOUBLE, i, 0, row_comm);
        }

        for (int i = 0; i < clusterCols; i++) {
            temp_row[i].clear();
            temp_col[i].clear();
            temp_data[i].clear();
        }
        temp_row.clear();
        temp_col.clear();
        temp_data.clear();
    }

    // column masters recieve data from cluster master
    if (myId < clusterRows && myId != 0) {
        std::cout << myId << "column master reading data" << std::endl;
        MPI_Recv(&nodeColCount, 1, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);
        MPI_Recv(&colRowCount, 1, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);
        colMasterTemp_col.resize(nodeColCount);
        colMasterTemp_data.resize(nodeColCount);
        colMasterTemp_row.resize(colRowCount);

        MPI_Recv(&colMasterTemp_row[0], colRowCount, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);
        MPI_Recv(&colMasterTemp_col[0], nodeColCount, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);
        MPI_Recv(&colMasterTemp_data[0], nodeColCount, MPI_DOUBLE, 0, 0, row_comm, MPI_STATUS_IGNORE);
    }

    if (myId < clusterRows) {
        for (int k = 0; k < rowsPerNode; k++) {
            csr_row.push_back(csr_data.size());

            for (int l = colMasterTemp_row[k]; l < colMasterTemp_row[k+1]; l++) {
                csr_data.push_back(colMasterTemp_data[l]);
                csr_col.push_back(colMasterTemp_col[l]);
            }
        }
    }

    // column masters to send column nodes their individual data
    if (myId < clusterRows) {
        for (int i = 1; i < clusterRows; i++){
            if (colMasterTemp_data.size() == 0){
                int nothingToSend = 0;
                MPI_Send(&nothingToSend, 1, MPI_INT, i, 0, col_comm);
            } else {
                int j = 0;
                while (colMasterTemp_row[(i * rowsPerNode) + j] == -1) { ;
                    j++;
                }

                int firstDataAt = colMasterTemp_row[(i * rowsPerNode) +j];
                int lastDataAt;
                if (i == clusterRows - 1) {
                    lastDataAt = colMasterTemp_col.size();
                } else {
                    int k = 0;
                    while (colMasterTemp_row[((i + 1) * (rowCount / clusterRows)) - k] == -1) {
                        k++;
                    }
                    lastDataAt = colMasterTemp_row[((i + 1) * (rowCount / clusterRows)) - k];
                }

                int numElements = (lastDataAt - firstDataAt);
                MPI_Send(&numElements, 1, MPI_INT, i, 0, col_comm);

                if (numElements > 0) {
                    MPI_Send(&colMasterTemp_row[i * (rowCount / clusterRows)], rowCount / clusterRows, MPI_INT, i, 0,
                             col_comm);
                    MPI_Send(&colMasterTemp_col[firstDataAt], numElements, MPI_INT, i, 0, col_comm);
                    MPI_Send(&colMasterTemp_data[firstDataAt], numElements, MPI_DOUBLE, i, 0, col_comm);
                }
            }
        }


    }

    if(myId >= clusterRows) {
        std::cout << "Non-master " << myId << " about to read data sent from column master" << std::endl;
        MPI_Recv(&nodeColCount, 1, MPI_INT, 0, 0, col_comm, MPI_STATUS_IGNORE);

            csr_col.resize(nodeColCount);
            csr_data.resize(nodeColCount);
            csr_row.resize(rowCount/clusterRows);

        if (nodeColCount > 0) {
            MPI_Recv(&csr_row[0], rowCount/clusterRows, MPI_INT, 0, 0, col_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&csr_col[0], nodeColCount, MPI_INT, 0, 0, col_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&csr_data[0], nodeColCount, MPI_DOUBLE, 0, 0, col_comm, MPI_STATUS_IGNORE);
        }
        std::cout << myId << "  DONE RECIEVING DATA " << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (!myId) std::cout << "Broadcasting dense vector to all nodes" << std::endl;
    MPI_Bcast(&denseVector[0], rowCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << myId << " About to calculate node results" << std::endl;

    double nodeResult[rowsPerNode] = {0.0};
    if (!myId){
        nodeColCount = csr_data.size();
    }

    if (nodeColCount > 0) { // check if this node has any data to work on at all
        clusterSpMV(ompThreads, csr_row, csr_col, csr_data, denseVector, nodeResult, rowsPerNode);

        csr_row.clear();
        csr_col.clear();
        csr_data.clear();
    }

    std::cout << myId << " - DONE CALCULATING" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    double *rowResult;
    if (((myId % clusterRows) == 0) || myId == 0) {
        rowResult = new double[rowsPerNode];
    }

    MPI_Reduce(&nodeResult, rowResult, rowsPerNode, MPI_DOUBLE, MPI_SUM, 0, row_comm);
    if ((myId % clusterRows) == 0) {
        if (!myId) {
            SpMVResult = new double[rowCount];
            std::cout << myId << "Master will gather " << rowsPerNode << " elements from row masters" << std::endl;
        }

        MPI_Barrier(col_comm);
        MPI_Gather(&rowResult[0], rowsPerNode, MPI_DOUBLE, &SpMVResult[0], elementsPerRow, MPI_DOUBLE, 0, col_comm);
    }

    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&row_comm);
    MPI_Finalize();

    if (!myId) {
        int miscalculations = 0;
        int localZeros = 0, distributedZeros = 0;
        for (int i = 0; i < rowCount; i++){
            if (SpMVResult[i] == 0.0) {
                distributedZeros++;
            }
            if (masterOnly_SpMV[i] == 0.0) {
                localZeros++;
            }
            if (SpMVResult[i] != masterOnly_SpMV[i]) {
                std::cout << "row " << i << ": " << SpMVResult[i] << " != " << masterOnly_SpMV[i] << std::endl;
                miscalculations++;
            }
            //std::cout << "row " << i << ": " << SpMVResult[i] << " and  " << masterOnly_SpMV[i] << std::endl;
        }

        if (localZeros) std::cout << "*** " << localZeros << " Local Zero Value Rows ***" << std::endl;
        if (distributedZeros) std::cout << "*** " << distributedZeros << " Distributed Zero Value Rows ***" << std::endl;
        if (miscalculations) std::cout << "*** " << miscalculations << " Miscalculations ***" << std::endl;
    }

    return 0;
}