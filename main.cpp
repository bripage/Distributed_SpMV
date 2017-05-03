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
    std::vector<int> origin_row;
    std::vector<int> origin_col;
    std::vector<double> origin_data;
    int rowCount, colCount, nonZeros;

    if (!myId) {
        //convert sparse matrix from Matrix Market to Compressed Sparse Row format
        MMCOO_to_CSR(argv[1], origin_row, origin_col, origin_data, rowCount, colCount, nonZeros);
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
                    std::cout << "temp += origin_data[" << j << "] * denseVector[origin_col[" << j << "]] --> " << origin_data[j] << " * " << denseVector[origin_col[j]] << std::endl;
                }
            } else {
                for (int j = origin_row[i]; j < origin_data.size(); j++) {
                    temp += origin_data[j] * denseVector[origin_col[j]];
                    std::cout << "temp += origin_data[" << j << "] * denseVector[origin_col[" << j << "]] --> " << origin_data[j] << " * " << denseVector[origin_col[j]] << std::endl;
                }
            }
            masterOnly_SpMV[i] = temp;
            std::cout << "masterOnly_SpMV[" << i << "] = " << temp << std::endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // create arrays for individual nodes to hold their data.
    int nodeColCount, colRowCount;
    int rowsPerNode = rowCount / clusterRows;
    int colsPerNode = rowsPerNode;

    std::vector<int> csr_row, csr_col, colMasterTemp_row, colMasterTemp_col;
    std::vector<double> csr_data, colMasterTemp_data;


    if (!myId) {
        std::vector<std::vector<int> > temp_row,temp_col;
        std::vector<std::vector<double> > temp_data;
        for (int i = 0; i < clusterCols; i++) {
            std::vector<int> a,b;
            std::vector<double> c;
            temp_row.push_back(a);
            temp_col.push_back(b);
            temp_data.push_back(c);
        }
//std::cout << "1" << std::endl;
        /*
        for (int i = 0; i < rowCount; i++) {
            int rowStart = origin_row[i];
            int nextRowAt;
            std::cout << "2" << std::endl;
            if (i >= rowCount) {
                nextRowAt = -1;
            } else {
                nextRowAt = origin_row[i+1];
            }
            std::cout << "3" << std::endl;
            int elementCount;
            if (nextRowAt == -1) {
                elementCount = origin_data.size() - rowStart;
            } else {
                elementCount = nextRowAt - rowStart;
            }
            std::cout << "4" << std::endl;
            if (elementCount < colsPerNode) {
                std::cout << "inside if" << std::endl;
                std::cout << "temp_row[0].push_back(" << temp_col[0].size() << ")" << std::endl;

                temp_row[0].push_back(temp_col[0].size());
                std::cout << "5" << std::endl;
                for (int k = 1; k < clusterRows; k++) {
                    temp_row[k].push_back(-1);
                }
                std::cout << "6" << std::endl;
            } else {
                std::cout << "7" << std::endl;
                std::cout << "inside else" << std::endl;
                int fullClusterCols = elementCount / colsPerNode;
                std::cout << "fullRows = " << fullClusterCols << std::endl;
                for (int k = 0; k < fullClusterCols+1; k++) {
                    std::cout << "8" << std::endl;
                    std::cout << "inside this loop" << std::endl;
                    temp_row[k].push_back(temp_col[k].size());
                    if(k == 0){
                        std::cout << "temp_row[0].push_back(" << temp_col[k].size() << ")" << std::endl;
                    }
                    std::cout << "9" << std::endl;
                }
                for (int k = fullClusterCols+1; k < rowsPerNode; k++) {
                    temp_row[k].push_back(-1);
                    std::cout << "10" << std::endl;
                }
            }
            std::cout << "11" << std::endl;
            if (nextRowAt != -1){
                std::cout << "12" << std::endl;
                //go until the end of the current row and no further in this iteration
                for (int j = rowStart; j < nextRowAt; j++){
                    std::cout << "13" << std::endl;
                    std::cout << "temp_col[(" << i << "/" << colsPerNode << ")].push_back(origin_col[" << j
                              << "]) --> temp_col[" << (i/colsPerNode) << "].push_back(" << origin_col[j]
                              << ")" << std::endl;
                        temp_col[(i/colsPerNode)].push_back(origin_col[j]);

                    std::cout << "temp_data[(" << i << "/" << colsPerNode << ")].push_back(origin_data[" << j
                              << "]) --> temp_data[" << (i/colsPerNode) << "].push_back(" << origin_data[j]
                              << ")" << std::endl;
                        temp_data[(i/colsPerNode)].push_back(origin_data[j]);
                }
            } else {
                std::cout << "15" << std::endl;
                // go till end of temp_col.end()
                for (int j = rowStart; j < origin_col.size(); j++){
                    std::cout << "16" << std::endl;
                    temp_col[j/colsPerNode].push_back(origin_col[j]);
                    temp_data[j/colsPerNode].push_back(origin_data[j]);
                    std::cout << "17" << std::endl;
                }
            }
        }

        for (int k = 0; k < temp_data[0].size(); k++){
            csr_data.push_back(temp_data[0][k]);
            csr_col.push_back(temp_col[0][k]);
            std::cout << "18" << std::endl;
        }
        for (int k = 0; k < rowsPerNode; k++) {
            csr_row.push_back(temp_row[0][k]);
            std::cout << "19" << std::endl;
        }

        for (int k = 0; k < temp_row.size(); k++){
            //std::cout << temp_row[0][k] << std::endl;
        }
*/
        for (int k = 0; k < origin_row.size(); k++){
            std::cout << origin_row[k] << std::endl;
        }

        for (int i = 0; i < rowCount; i++) {
            int firstData = origin_row[i];
            int lastData, rowElements;

            if (i == rowCount - 1){
                lastData = origin_data.size()-1;
                rowElements = origin_data.size() - firstData;
            } else {
                lastData = origin_row[i+1]-1;
                rowElements = lastData - firstData;
            }

            std::cout << "i = " << i << ", firstData = " << firstData << ", lastData = " << lastData << std::endl;
            int maxMatCols = colsPerNode;

            int clusterColsNeeded = rowElements / maxMatCols;
            std::cout << "clusterColsNeeded = " << clusterColsNeeded << std::endl;
            for (int k = 0; k < clusterCols; k++){
                if (k <= clusterColsNeeded) {
                    temp_row[k].push_back(temp_data[k].size());
                    std::cout << "temp_row[" << k << "].push_back(" << temp_data[k].size() << ")" << std::endl;
                } else {
                    temp_row[k].push_back(-1);
                    std::cout << "temp_row[" << k << "].push_back(-1)" << std::endl;
                }
            }

            int count = 0;
            for (int j = firstData; j <= lastData; j++){
                temp_data[count / maxMatCols].push_back(origin_data[j]);
                //std::cout << "temp_data[" << count << "/" <<  maxMatCols << "].push_back(origin_data[" << j << "]) ---> "
                //          << "temp_data[" << count / maxMatCols << "].push_back(" << origin_data[j] << ")" << std::endl;
                temp_col[count / maxMatCols].push_back(origin_col[j]);

                //std::cout << "count = " << count << ", csrDataSize = " << csr_data.size() << ", tempDataSize = " << temp_data[0].size() << std::endl;
                count++;
            }
        }

        for (int k = 0; k < temp_data[0].size(); k++){
            colMasterTemp_data.push_back(temp_data[0][k]);
            colMasterTemp_col.push_back(temp_col[0][k]);
            //std::cout << "18" << std::endl;
        }
        for (int k = 0; k < rowCount; k++) {
            colMasterTemp_row.push_back(temp_row[0][k]);
            //std::cout << "19" << std::endl;
        }

       //std::cout << "Before clearing temps -> " << temp_row.size() << ", " << temp_col.size() << ", " << temp_data.size() << std::endl;

        //MPI_Barrier(MPI_COMM_WORLD);

        // Send column data to column masters
        for ( int i = 1; i < clusterRows; i++) {
                nodeColCount = temp_data[i].size();
                colRowCount = temp_row[i].size();
                //std::cout << temp_row.size() << ", " << temp_col.size() << ", " << temp_data.size() << std::endl;
                //std::cout << temp_row[i].size() << ", " << temp_col[i].size() << ", " << temp_data[i].size() << std::endl;
                //std::cout << rowCount << std::endl;

                //std::cout << "Master sending nodeColCount to " << i << std::endl;
                MPI_Send(&nodeColCount, 1, MPI_INT, i, 0, row_comm);
                //MPI_Send(&rowCount, 1, MPI_INT, i, 0, row_comm);
                MPI_Send(&colRowCount, 1, MPI_INT, i, 0, row_comm);
                //std::cout << "Master sending rows to " << i << std::endl;
                MPI_Send(&temp_row[i][0], rowCount, MPI_INT, i, 0, row_comm);
                //std::cout << "Master sending cols " << i << std::endl;
                MPI_Send(&temp_col[i][0], nodeColCount, MPI_INT, i, 0, row_comm);
                //std::cout << "Master sending data to " << i << std::endl;
                MPI_Send(&temp_data[i][0], nodeColCount, MPI_DOUBLE, i, 0, row_comm);

                // Free up column data
                //temp_row[i].clear();
                //temp_col[i].clear();
                //temp_data[i].clear();
        }
        //std::cout << "After clearing temps -> " << temp_row.size() << ", " << temp_col.size() << ", " << temp_data.size() << std::endl;


        //std::cout << "Master about to clear CSR temp vectors" << std::endl;
        origin_row.clear();
        origin_col.clear();
        origin_data.clear();
        //std::cout << "Master done clearing CSR temp vectors" << std::endl;
    }



    // column masters recieve data from cluster master
    if (myId < clusterRows && myId != 0) {
        std::cout << myId << "column master reading data" << std::endl;
        MPI_Recv(&nodeColCount, 1, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);
        MPI_Recv(&colRowCount, 1, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);
        //MPI_Recv(&rowCount, 1, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);
        //std::cout << "Non-master recieved col count" << std::endl;
        colMasterTemp_col.resize(nodeColCount);
        colMasterTemp_data.resize(nodeColCount);
        colMasterTemp_row.resize(colRowCount);
        //std::cout <<  myId << "rowCount = " << rowCount << std::endl;
        //std::cout << myId << " resized vectors" << std::endl;
        //std::cout << myId <<": " << csr_row.size() << ", " << csr_col.size() << ", " << csr_data.size() << std::endl;

        //std::cout << myId << "  about to recieve rows" << std::endl;
        MPI_Recv(&colMasterTemp_row[0], colRowCount, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);
        //std::cout << myId << "  recieved rows" << std::endl;
        MPI_Recv(&colMasterTemp_col[0], nodeColCount, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);
        //std::cout << myId << "  recieved cols" << std::endl;
        MPI_Recv(&colMasterTemp_data[0], nodeColCount, MPI_DOUBLE, 0, 0, row_comm, MPI_STATUS_IGNORE);
        //std::cout << myId << "  received data" << std::endl;
    }

    if (myId < clusterRows) {
        for (int k = 0; k < rowsPerNode; k++) {
            csr_row.push_back(csr_data.size());
            std::cout << myId << "- csr_row[" << k << "] = " << csr_row[k] << std::endl;

            for (int l = colMasterTemp_row[k]; l < colMasterTemp_row[k+1]; l++) {
                std::cout << myId << "- csr_data.push_back(" << colMasterTemp_data[l] << ")" << std::endl;
                csr_data.push_back(colMasterTemp_data[l]);
                csr_col.push_back(colMasterTemp_col[l]);
                //std::cout << "18" << std::endl;
            }
        }
    }

    if (myId < clusterRows) {
        std::cout << myId << " - colMasterTemp_data.size() = " << colMasterTemp_data.size() << std::endl;
    }
    //MPI_Barrier(MPI_COMM_WORLD);



    // column masters to send column nodes their individual data
    if (myId < clusterRows) {
        for (int i = 1; i < clusterRows; i++){
            if (colMasterTemp_data.size() == 0){
                int nothingToSend = 0;
                MPI_Send(&nothingToSend, 1, MPI_INT, i, 0, col_comm);
            } else {
                int j = 0;
                std::cout << myId << ": " << colMasterTemp_row.size() << " rows, " << colMasterTemp_col.size()
                          << " cols, "
                          << colMasterTemp_data.size() << " data" << std::endl;
                while (colMasterTemp_row[(i * rowsPerNode) + j] == -1) {
                    std::cout << i << " * " << (rowCount / clusterRows) << " + " << j << " = "
                              << (i * (rowCount / clusterRows)) + j << ", "
                              << colMasterTemp_row[(i * (rowCount / clusterRows)) + j]
                              << std::endl;
                    j++;
                }
                int firstDataAt = colMasterTemp_row[(i * rowsPerNode) +j];
                std::cout << i * clusterRows << " firstDataAt = " << firstDataAt << std::endl;
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

                std::cout << i * clusterRows << " lastDataAt = " << lastDataAt << std::endl;
                int numElements = (lastDataAt - firstDataAt);

                std::cout << myId << ": sending numElements " << numElements << " to " << i << std::endl;
                MPI_Send(&numElements, 1, MPI_INT, i, 0, col_comm);
                //MPI_Send(&rowCount, 1, MPI_INT, i, 0, col_comm);

                if (numElements > 0) {
                    std::cout << "sendings rows to " << i << std::endl;
                    MPI_Send(&colMasterTemp_row[i * (rowCount / clusterRows)], rowCount / clusterRows, MPI_INT, i, 0,
                             col_comm);
                    std::cout << "sending cols to " << i << std::endl;
                    MPI_Send(&colMasterTemp_col[firstDataAt], numElements, MPI_INT, i, 0, col_comm);
                    std::cout << "sending data to " << i << std::endl;
                    MPI_Send(&colMasterTemp_data[firstDataAt], numElements, MPI_DOUBLE, i, 0, col_comm);
                }
            }
        }


    }

    if(myId >= clusterRows) {
        std::cout <<  myId << "Non-master about to read data sent from column master" << std::endl;
        MPI_Recv(&nodeColCount, 1, MPI_INT, 0, 0, col_comm, MPI_STATUS_IGNORE);
        //MPI_Recv(&rowCount, 1, MPI_INT, 0, 0, col_comm, MPI_STATUS_IGNORE);
        std::cout <<  myId << "Non-master recieved col count of " << nodeColCount << std::endl;

        //if (nodeColCount > 0) {
            csr_col.resize(nodeColCount);
            csr_data.resize(nodeColCount);
            csr_row.resize(rowCount/clusterRows);

            std::cout << myId << " resized vectors" << std::endl;
        if (nodeColCount > 0) {
            MPI_Recv(&csr_row[0], rowCount/clusterRows, MPI_INT, 0, 0, col_comm, MPI_STATUS_IGNORE);
            std::cout << myId << "  recieved rows" << std::endl;
            MPI_Recv(&csr_col[0], nodeColCount, MPI_INT, 0, 0, col_comm, MPI_STATUS_IGNORE);
            std::cout << myId << "  recieved cols" << std::endl;
            MPI_Recv(&csr_data[0], nodeColCount, MPI_DOUBLE, 0, 0, col_comm, MPI_STATUS_IGNORE);
            std::cout << myId << "  received data" << std::endl;
        }
        std::cout << myId << "  DONE RECIEVING DATA " << std::endl;
    }


    if (myId == 0 || myId == 2){
        for (int i = 0; i < csr_row.size(); i ++){
            std::cout << myId << "- csr_row[" << i << "] = " << csr_row[i] << std::endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (!myId) std::cout << "Broadcasting dense vector to all nodes" << std::endl;
    //std::cout << "Broadcasting dense vector to all nodes" << std::endl;
    MPI_Bcast(&denseVector[0], rowCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << myId << " About to calculate node results" << std::endl;
    if (nodeColCount > 0) {
        std::cout << myId << "Has " << nodeColCount << " cols to work on" << std::endl;
    }

    double nodeResult[rowsPerNode] = {0.0};
    std::cout << myId << " Done creating nodeResult array" << std::endl;

    if (!myId){
        nodeColCount = csr_data.size();
    }
    if (nodeColCount > 0) { // check if this node has any data to work on at all
        std::cout << myId << " is inside if 1" << std::endl;
        //for (int i = 0; i < csr_row.size(); i++) { // iterate through all rows node is to work on
        for (int i = 0; i < rowsPerNode; i++) { // iterate through all rows node is to work on
            std::cout << myId << " is inside for 1" << std::endl;
            if (csr_row[i] != -1) { // check is i points to row that has no data on this node
                std::cout << myId << " is inside if 2" << std::endl;
                if (i == (rowsPerNode) - 1) { // check if i is last row of node
                    std::cout << myId << " is inside if 3" << std::endl;
                    for (int j = csr_row[i]-csr_row[0]; j < csr_col.size(); j++) { // go from last row start to end of data
                        std::cout << myId << " is inside for 2" << std::endl;
                        nodeResult[i] += csr_data[j] * denseVector[csr_col[j]];
                        std::cout << myId << "nodeResult[" << i << "] += csr_data[" << j << "] * denseVector[csr_col[" << j << "]] --> " << csr_data[j] << " * " << denseVector[csr_col[j]] << std::endl;
                    }
                }else {
                    int nextValidRow = i+1;
                    while ((csr_row[nextValidRow] == -1) && (nextValidRow < csr_row.size())){
                        nextValidRow++;
                    }
                    std::cout << "nextValidRow = " << nextValidRow << std::endl;

                    int rowStart = csr_row[i]-csr_row[0];
                    int rowEnd = csr_row[nextValidRow]-csr_row[0];

                    for (int j = rowStart; j < rowEnd; j++) { // go from last row start to before row i+1
                        nodeResult[i] += csr_data[j] * denseVector[csr_col[j]];
                        std::cout << myId << "nodeResult[" << i << "] += csr_data[" << j << "] * denseVector[csr_col[" << j << "]] --> " << csr_data[j] << " * " << denseVector[csr_col[j]] << std::endl;
                    }
                }
            }
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << myId << " - DONE CALCULATING" << std::endl;

    double *rowResult;
    if (((myId % clusterRows) == 0) || myId == 0) {
        rowResult = new double[elementsPerRow];
        std::cout << myId << " created rowResult array of size " << elementsPerRow << std::endl;
    }
    std::cout << myId << " 20" << std::endl;
        //MPI_Barrier(col_comm);
        MPI_Reduce(&nodeResult, rowResult, elementsPerRow, MPI_DOUBLE, MPI_SUM, 0, row_comm);
    std::cout << myId << " 21" << std::endl;
    if ((myId % clusterRows) == 0) {
        if (!myId) {
            SpMVResult = new double[numToGather];
            std::cout << myId << " to gather " << numToGather << "elements from row masters" << std::endl;
        }

        //MPI_Barrier(col_comm);
        MPI_Gather(&rowResult[0], elementsPerRow, MPI_DOUBLE, &SpMVResult[0], elementsPerRow, MPI_DOUBLE, 0, col_comm);
    }

    std::cout << myId << " 22" << std::endl;

    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&row_comm);
    MPI_Finalize();

    if (!myId) {
        std::cout << "SpMV Result contains " << numToGather << " values" << std::endl;

        if (numToGather != rowCount) {
            std::cout << "SOME ROWS LOST DURING DISTRIBUTION (Non-Zeros Integer Division Problem)" << std::endl;
        }

        int miscalculations = 0;
        int localZeros = 0, distributedZeros = 0;
        for (int i = 0; i < numToGather; i++){
            if (SpMVResult[i] == 0.0) {
                distributedZeros++;
            }
            if (masterOnly_SpMV[i] == 0.0) {
                localZeros++;
            }
            if (SpMVResult[i] != masterOnly_SpMV[i]) {
                miscalculations++;
            }
            std::cout << i << ": " << SpMVResult[i] << " != " << masterOnly_SpMV[i] << std::endl;
        }
        if (localZeros) std::cout << "*** " << localZeros << " Local Zero Value Rows ***" << std::endl;
        if (distributedZeros) std::cout << "*** " << distributedZeros << " Distributed Zero Value Rows ***" << std::endl;
        if (miscalculations) std::cout << "*** " << miscalculations << " Miscalculations ***" << std::endl;
    }


    return 0;
}