#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <iomanip>
#include "csrClusterSplit.h"
#include "clusterSpMV.h"
#include "mm_to_csr.h"
#include <string>

int main(int argc, char *argv[]) {

    std::string argTemp, distributionMethod;
    bool masterOnly = false, colMajor = false;
    int ompThreads = 1;
    char *matrixFile, *vectorFile;
    for (int i = 1; i < argc; i= i+2){
        argTemp = argv[i];
        if (argTemp == "-lm"){
            // load matrix from file.
            matrixFile = argv[i+1];
        } else if (argTemp == "-lv"){
            // load dense vector from file.
            vectorFile = argv[i+1];
        } else if (argTemp == "-s"){
            //run on MPI master only (no distribution, OpenMP as normal).
            if (argv[i+1] == "true"){
                masterOnly = true;
            }
        } else if (argTemp == "--distribution-method"){
            //set number of OpenMP threads per node
            distributionMethod = argv[i+1];
        } else if (argTemp == "--major-order") {
            std::string temp = argv[i + 1];
            if (temp == "col") {
                colMajor = true;
            }
        } else if (argTemp == "--omp-threads"){
                //set number of OpenMP threads per node
                ompThreads = atoi(argv[i+1]);
        } else if (argTemp == "--help"){
            std::cout << "Usage: distSpMV [OPTION] <argument> ..." << std::endl << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << " -lm <file>" << std::setw(15) << "" << R"(Load sparse matrix from <file>)" << std::endl;
            std::cout << " -lv <file>" << std::setw(15) << "" << R"(Load dense vector from <file>)" << std::endl;
            std::cout << " -s <arg>" << std::setw(15) << "" << R"(Set master only computation. In master only
            computation no distribution to cluster members will occur, however OpenMP will still operate with the set
            number of threads per node)" << std::endl;

            std::cout << " --distribution-method <arg>" << std::setw(15) << "" << R"(Set the work distribution
            algorotithm to be used for partitioning work amongst cluster nodes)" << std::endl;
            std::cout << " --major-order <arg>" << std::setw(15) << "" << R"(Set the major order of the input matrix, as
            well as compression and calculations. This application supports both row and col major order. )"
                      << std::endl;
            std::cout << " --omp-threads <arg>" << std::setw(15) << "" << R"(Set the number of OpenMP threads per node)"
                      << std::endl;

            exit(0);
        } else {
            printf("%s Is not a valid parameter. EXITING!\n", argv[i]);
            exit(0);
        }
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

    // Split dupCommWorld into comm's for each node column, based on each processes original myIdhttps://duckduckgo.com/?t=canonical&atb=v53-4__
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
    std::vector<int> origin_row, csr_row, origin_col, csr_col;
    std::vector<double> origin_data, csr_data;
    std::vector<std::vector<int> > temp_row, temp_col,colMasterTemp_row, colMasterTemp_col;
    std::vector<std::vector<double> > temp_data, colMasterTemp_data;

    double *masterOnly_SpMV, temp;
    std::vector<double> denseVector;

    int rowCount, colCount, nonZeros;
    int colsPerNode;

    if (!myId) {
        //convert sparse matrix from Matrix Market to Compressed Sparse Row format
        if (colMajor) {
            std::cout << "*** Using column major order of input matrix ***" << std::endl;
            MMCOO_to_CSR_colMajor(matrixFile, origin_row, origin_col, origin_data, rowCount, colCount, nonZeros);
        } else {
            MMCOO_to_CSR(matrixFile, origin_row, origin_col, origin_data, rowCount, colCount, nonZeros);
        }

        //std::cout << "Creating Dense Vector" << std::endl;
        if (!myId) {
            //std::cout << "Populating Dense Vector" << std::endl;
            for (int i = 0; i < rowCount; i++) {
                denseVector.push_back(1.0);
            }
        }

        std::cout << "Performing Master Only SpMV" << std::endl;
        masterOnly_SpMV = new double[rowCount];
        for (int i = 0; i < rowCount; i++) {
            temp = 0.0;
            if (i != rowCount - 1) {
                if (colMajor) { //col major order selected
                    for (int j = origin_row[i]; j < origin_row[i + 1]; j++) {   // go to end of current row
                        // entire row is multiplied by a single dense vector element
                        temp += origin_data[j] * denseVector[i];
                    }
                } else {    // row major order selected
                    for (int j = origin_row[i]; j < origin_row[i + 1]; j++) {   // go to end of current row
                        temp += origin_data[j] * denseVector[origin_col[j]];
                    }
                }
            } else {
                if (colMajor) { //col major order selected
                    for (int j = origin_row[i]; j < origin_data.size(); j++) {  // go to end of data vector
                        // entire row is multiplied by a single dense vector element
                        temp += origin_data[j] * denseVector[i];
                    }
                } else {    // row major order selected
                    for (int j = origin_row[i]; j < origin_data.size(); j++) {  // go to end of data vector
                        temp += origin_data[j] * denseVector[origin_col[j]];
                    }
                }
            }
            masterOnly_SpMV[i] = temp;
        }
    }


    //***********************************************//
    //     Select distribution method and actions    //
    //***********************************************//
    MPI_Barrier(MPI_COMM_WORLD);

    if (distributionMethod == "splitmatrix") {
        std::vector<std::vector <std::vector <int> > > nodeRowOwnership;
        std::vector<std::vector <double> > splitDenseVector;
        int colsLastColumn, colsPerColumn, nodeElementCount, nodeRowCount, colElementCount;

        if (!myId) {
            colsPerColumn = rowCount / clusterCols;
            if (rowCount % clusterCols != 0) {
                colsLastColumn = colsPerColumn + (rowCount % clusterCols);
            } else {
                colsLastColumn = rowCount / clusterCols;
            }

            csrClusterSplit_SplitMatrix(matrixFile, colMajor, distributionMethod, origin_row, origin_col, origin_data,
                                        colMasterTemp_row, colMasterTemp_col, colMasterTemp_data, rowCount, colCount,
                                        nonZeros, colsPerNode, clusterRows, clusterCols);


            std::cout << "Distribution of non-zero elements:" << std::endl;
            //std::cout << "temp.size() = " << temp_data.size() << std::endl;
            for (int i = 0; i < clusterCols; i++) {
                std::cout << "Column " << i << ": " << colMasterTemp_data[i].size() << " elements assigned"
                          << std::endl;
                for (int j = 0; j < clusterRows; j++) {
                    if (colMasterTemp_data[i].size() == 0) {
                        std::cout << "\tRow " << j << ": 0" << std::endl;
                    } else {
                        int start, stop, k = 0, l = 0;
                        while (colMasterTemp_row[i][(j * colsPerColumn) + k] == -1) {
                            k++;
                        }
                        start = colMasterTemp_row[i][(j * colsPerColumn) + k];

                        if (j == clusterRows - 1) {
                            stop = colMasterTemp_data[i].size();
                        } else {
                            while (colMasterTemp_row[i][((j + 1) * colsPerColumn) + l] == -1) {
                                l++;
                            }
                            stop = colMasterTemp_row[i][((j + 1) * colsPerColumn) + l];
                        }

                        std::cout << "\tRow " << j << ": " << stop - start << std::endl;
                    }
                }
            }
            std::cout << std::endl;

            // send column master data to individual column masters
            for (int i = 1; i < clusterCols; i++) {
                colElementCount = colMasterTemp_data[i].size(); // how many elements clusterCol i has to work with
                std::cout << "colElementCount = " <<  colElementCount << std::endl;
                std::cout << "colMasterTemp_row[" << i << "].size() = " <<  colMasterTemp_row[i].size() << ", rowCount = " << rowCount << std::endl;

                MPI_Send(&colElementCount, 1, MPI_INT, i, 0, row_comm);
                MPI_Send(&rowCount, 1, MPI_INT, i, 0, row_comm);
                MPI_Send(&colMasterTemp_row[i][0], rowCount, MPI_INT, i, 0, row_comm);
                //MPI_Send(&colMasterTemp_col[i][0], colElementCount, MPI_INT, i, 0, row_comm);
                //MPI_Send(&colMasterTemp_data[i][0], colElementCount, MPI_DOUBLE, i, 0, row_comm);

                // Dont need to keep this data on the column master once its sent to the proper node!
                //colMasterTemp_row[i].clear();
                //colMasterTemp_col[i].clear();
                //colMasterTemp_data[i].clear();
            }
        }

        // column masters recieve data from cluster master
        if (myId < clusterRows && myId != 0) {
            std::cout << myId << "column master reading data" << std::endl;
            MPI_Recv(&colElementCount, 1, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);  // get number of elements
            MPI_Recv(&rowCount, 1, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE); // get number of rows (should be all)

            std::cout << "rowCount = " << rowCount << std::endl;
            // resize the vectors to allow for the incoming data to be recieved
            colMasterTemp_col.resize(colElementCount);
            colMasterTemp_data.resize(colElementCount);
            colMasterTemp_row.resize(rowCount);

            MPI_Recv(&colMasterTemp_row[0], rowCount, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&colMasterTemp_col[0], colElementCount, MPI_INT, 0, 0, row_comm, MPI_STATUS_IGNORE);
           //MPI_Recv(&colMasterTemp_data[0], colElementCount, MPI_DOUBLE, 0, 0, row_comm, MPI_STATUS_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Comm_free(&col_comm);
        MPI_Comm_free(&row_comm);
        MPI_Finalize();

    } else if (distributionMethod == "colbalanced"){
        std::vector<std::vector <std::vector <int> > > nodeRowOwnership;
        std::vector<std::vector <double> > splitDenseVector;

        if (!myId) {
            csrClusterSplit_ElementBalanced(matrixFile, colMajor, distributionMethod, processCount, origin_row,
                                            origin_col, origin_data, temp_row, temp_col, temp_data, rowCount, colCount,
                                            nonZeros, colsPerNode, clusterRows, clusterCols, nodeRowOwnership);

            splitDenseVector_ElementBalanced(denseVector, splitDenseVector, nodeRowOwnership);

            std::cout << "Distribution of non-zero elements:" << std::endl;
            for (int i = 0; i < processCount; i++){
                std::cout << "Node " << i << ": " << nodeRowOwnership[i][0][0] << std::endl;
            }
        }


        // perform SpMV on master only, without any segmentation or distribution
        // for verification purposes.
        if (!myId) {
            origin_row.clear();
            origin_col.clear();
            origin_data.clear();
        }


        int nodeElementCount, nodeRowCount;
        int rowsPerNode = rowCount / clusterRows;
        //std::vector<int> myRows;

        if (!myId) {
            // give master its work data
            csr_row = temp_row[0];
            csr_col = temp_col[0];
            csr_data = temp_data[0];

            // Send column data to nodes
            for (int i = 1; i < processCount; i++) {
                nodeElementCount = temp_data[i].size();
                nodeRowCount = temp_row[i].size();

                MPI_Send(&nodeRowCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&nodeElementCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                //MPI_Send(&nodeRowOwnership[i][2], nodeRowCount, MPI_INT, i, 0, row_comm);
                MPI_Send(&temp_row[i][0], nodeRowCount, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&temp_col[i][0], nodeElementCount, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&temp_data[i][0], nodeElementCount, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                MPI_Send(&splitDenseVector[i][0], nodeRowCount, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }

            for (int i = 0; i < clusterCols; i++) {
                temp_row[i].clear();
                temp_col[i].clear();
                temp_data[i].clear();
                splitDenseVector[i].clear();
            }
            temp_row.clear();
            temp_col.clear();
            temp_data.clear();
            splitDenseVector.clear();
        }

        /*
        if(!myId){
            // out put converted format
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << "----- Process " << myId << " has the following rows/data -----" << std::endl;
            for (int i = 0; i < nodeRowCount; i++){
                std::cout << "Row " << nodeRowOwnership[0][1][i] << ": ";
                if (i != nodeRowCount-1){
                    for (int j = csr_row[i]; j < csr_row[i+1]; j++){
                        std::cout << csr_data[j] << ", ";
                    }
                } else {
                    for (int j = csr_row[i]; j < csr_data.size(); j++){
                        std::cout << csr_data[j] << ", ";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
        */

        if (myId) {
            //std::cout << "Non-master " << myId << " about to read data " << std::endl;
            //std::cout << "0" << std::endl;
            MPI_Recv(&nodeRowCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //std::cout << "1" << std::endl;
            MPI_Recv(&nodeElementCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //std::cout << "2" << std::endl;
            //myRows.resize(nodeRowCount);
            csr_col.resize(nodeElementCount);
            csr_data.resize(nodeElementCount);
            csr_row.resize(nodeRowCount);
            denseVector.resize(nodeRowCount);
            //std::cout << "3" << std::endl;

            if (nodeElementCount > 0) {
                //MPI_Recv(&myRows[0], nodeRowCount, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&csr_row[0], nodeRowCount, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&csr_col[0], nodeElementCount, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&csr_data[0], nodeElementCount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&denseVector[0], nodeRowCount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            //std::cout << myId << "  DONE RECIEVING DATA " << std::endl;
        }

        std::vector <double> nodeResult;
        nodeResult.resize(nodeRowCount, 0.0);

        if (nodeElementCount > 0) { // check if this node has any data to work on at all
            //if (!myId) {
                //std::cout << "Calculating on " << myId << std::endl;
                clusterSpMV_ElementBalanced(ompThreads, csr_row, csr_col, csr_data, denseVector, nodeResult, colMajor);

                csr_row.clear();
                csr_col.clear();
                csr_data.clear();
            //}
        }
        //std::cout << "Done Calculating on " << myId << std::endl;

        // Send results to row master
        std::vector<std::vector <double> > clusterResults;
        //  Send back result to master
        if (!myId) {  // I am the master
            //copy master's local results into clusterResults
            clusterResults.push_back(nodeResult);

            std::cout << "Read results from " << myId << std::endl;
            for (int i = 1; i < processCount; i++){
                //std::cout << "Read results from " << i << std::endl;
                std::vector<double> tempResults(nodeRowOwnership[i][1].size());
                MPI_Recv(&tempResults[0], nodeRowOwnership[i][1].size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                clusterResults.push_back(tempResults);
            }
        } else {    // myId != 0, ie i am not the master
            //std::cout << "Send results from " << myId << std::endl;

            MPI_Send(&nodeResult[0], nodeRowCount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Comm_free(&col_comm);
        MPI_Comm_free(&row_comm);
        MPI_Finalize();

        if(!myId) {
            std::cout << "Done getting node data " << std::endl;

            double SpMVResult[rowCount] = {0.0};
            //  Combine results from slaves and any split rows
            for (int i = 0; i < processCount; i++) {
                for (int j = 0; j < nodeRowOwnership[i][1].size(); j++) {
                    //std::cout << "SpMVResult[" << nodeRowOwnership[i][1][j] << "] += clusterResults[" << i << "][" << j
                    //          << "]" << std::endl;
                    SpMVResult[nodeRowOwnership[i][1][j]] += clusterResults[i][j];
                }
            }

            std::cout << "Done combining results " << std::endl;

            int miscalculations = 0;
            int localZeros = 0, distributedZeros = 0;
            for (int i = 0; i < rowCount; i++) {
                if (SpMVResult[i] == 0.0) {
                    distributedZeros++;
                }
                if (masterOnly_SpMV[i] == 0.0) {
                    localZeros++;
                }
                if (SpMVResult[i] != masterOnly_SpMV[i] || (SpMVResult[i] == 0.0 || masterOnly_SpMV[i] == 0.0)) {
                    std::cout << "row " << i << ": " << SpMVResult[i] << " != " << masterOnly_SpMV[i] << std::endl;
                    miscalculations++;
                }
                //std::cout << "row " << i << ": " << SpMVResult[i] << " and  " << masterOnly_SpMV[i] << std::endl;
            }

            if (localZeros) std::cout << "*** " << localZeros << " Local Zero Value Rows ***" << std::endl;
            if (distributedZeros)
                std::cout << "*** " << distributedZeros << " Distributed Zero Value Rows ***" << std::endl;
            if (miscalculations) std::cout << "*** " << miscalculations << " Miscalculations ***" << std::endl;
        }
    }

    return 0;
}