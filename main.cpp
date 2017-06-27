#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include "controlStruct.h"
#include "csrSpMV.h"
#include "distribution.h"
#include "clusterSpMV.h"


int main(int argc, char *argv[]) {
    controlData controlData;
    std::string argTemp;

    for (int i = 1; i < argc; i= i+2){
        argTemp = argv[i];
        if (argTemp == "-lm"){
            // load matrix from file.
            controlData.matrixFile = argv[i+1];
        } else if (argTemp == "-lv"){
            // load dense vector from file.
            controlData.vectorFile = argv[i+1];
        } else if (argTemp == "-s"){
            //run on MPI master only (no distribution, OpenMP as normal).
            if (argv[i+1] == "true"){
                controlData.masterOnly = true;
            }
        } else if (argTemp == "--distribution-method"){
            //set number of OpenMP threads per node
            controlData.distributionMethod = argv[i+1];
        } else if (argTemp == "--major-order") {
            std::string temp = argv[i + 1];
            if (temp == "col") {
                controlData.colMajor = true;
            }
        } else if (argTemp == "--omp-threads"){
                //set number of OpenMP threads per node
            controlData.ompThreads = atoi(argv[i+1]);
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
    MPI_Comm_size(MPI_COMM_WORLD, &controlData.processCount);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &controlData.myId);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Verify MPI initialised properly
    if (controlData.processCount == 1 && controlData.masterOnly == false) {
        printf("MPI Initialization failed -- KILLING!");
        MPI_Finalize();
        exit(0);
    }

    //****************************************************//
    //  Create comm for each column/row of compute nodes  //
    //****************************************************//
    // get the number of MPI processes for work split and distribution
    controlData.clusterRows = sqrt(controlData.processCount);
    controlData.clusterCols = controlData.clusterRows; // works because cluster is expected to be square
    controlData.myCol = controlData.myId % controlData.clusterRows;

    // Create a duplicate of MPI_COMM_WORLD that can be used to split
    MPI_Comm dupCommWorldCol;
    MPI_Comm_dup(MPI_COMM_WORLD, &dupCommWorldCol);

    // Split dupCommWorld into comm's for each node column, based on each processes original myId
    //controlData.row_comm = new MPI_Comm;
    MPI_Comm_split(dupCommWorldCol, controlData.myCol, controlData.myId, &controlData.col_comm);

    //int col_rank, col_size;
    MPI_Comm_rank(controlData.col_comm, &controlData.col_rank);
    MPI_Comm_size(controlData.col_comm, &controlData.col_size);

    // get the number of MPI processes for work split and distribution
    controlData.myRow = controlData.myId / controlData.clusterRows;
    // Create a duplicate of MPI_COMM_WORLD that can be used to split
    MPI_Comm dupCommWorldRow;
    MPI_Comm_dup(MPI_COMM_WORLD, &dupCommWorldRow);

    // Split dupCommWorld into comm's for each node column, based on each processes original myId
    //controlData.row_comm = new MPI_Comm;
    MPI_Comm_split(dupCommWorldRow, controlData.myRow, controlData.myId, &controlData.row_comm);

    MPI_Comm_rank(controlData.row_comm, &controlData.row_rank);
    MPI_Comm_size(controlData.row_comm, &controlData.row_size);

    // create vectors to hold sparse matrix data once converted to CSR format
    std::vector<int> origin_row, csr_row, origin_col, csr_col;
    std::vector<double> origin_data, csr_data;
    std::vector<std::vector<int> > temp_row, temp_col,colMasterTemp_row, colMasterTemp_col;
    std::vector<std::vector<double> > temp_data, colMasterTemp_data;
    std::vector<double> denseVector;


    csrSpMV masterData;
    if (!controlData.myId) {
        masterData.masterOnlySpMV(controlData);
    }


    //***********************************************//
    //     Select distribution method and actions    //
    //***********************************************//
    MPI_Barrier(MPI_COMM_WORLD);

        std::vector<std::vector <std::vector <int> > > nodeRowOwnership;
        std::vector<std::vector <double> > splitDenseVector;
        std::vector<csrSpMV*> clusterColData;
        int colsPerColumn, colElementCount;

    if (!controlData.myId) {
        distribution_SplitMatrix(controlData, clusterColData);

    }


/*
        if (!controlData.myId) {
            colsPerColumn = controlData.rowCount / controlData.clusterCols;

            distribution_SplitMatrix(controlData, clusterColData);

            std::cout << "Distribution of non-zero elements:" << std::endl;
            //std::cout << "temp.size() = " << temp_data.size() << std::endl;
            for (int i = 0; i < controlData.clusterCols; i++) {
                std::cout << "Column " << i << ": " << colMasterTemp_data[i].size() << " elements assigned"
                          << std::endl;
                for (int j = 0; j < controlData.clusterRows; j++) {
                    if (colMasterTemp_data[i].size() == 0) {
                        std::cout << "\tRow " << j << ": 0" << std::endl;
                    } else {
                        int start, stop, k = 0, l = 0;
                        while (colMasterTemp_row[i][(j * colsPerColumn) + k] == -1) {
                            k++;
                        }
                        start = colMasterTemp_row[i][(j * colsPerColumn) + k];

                        if (j == controlData.clusterRows - 1) {
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

            origin_row.clear();
            origin_col.clear();
            origin_data.clear();


            // send column master data to individual column masters
            for (int i = 1; i < controlData.clusterCols; i++) {
                colElementCount = colMasterTemp_data[i].size(); // how many elements clusterCol i has to work with
                std::cout << "colElementCount = " <<  colElementCount << std::endl;
                std::cout << "colMasterTemp_row[" << i << "].size() = " <<  colMasterTemp_row[i].size() << ", rowCount = " << controlData.rowCount << "colMasterTemp_col[" << i << "].size() = " <<  colMasterTemp_col[i].size() << std::endl;

                MPI_Send(&colElementCount, 1, MPI_INT, i, 0, controlData.row_comm);
                std::cout << "Sent colElementCount" << std::endl;
                MPI_Send(&controlData.rowCount, 1, MPI_INT, i, 0, controlData.row_comm);
                std::cout << "Sent rowCount" << std::endl;
                MPI_Send(&(colMasterTemp_row[i][0]), controlData.rowCount, MPI_INT, i, 0, controlData.row_comm);
                std::cout << "Sent rows" << std::endl;
                MPI_Send(&colMasterTemp_col[i][0], colElementCount, MPI_INT, i, 0, controlData.row_comm);
                std::cout << "Sent cols" << std::endl;
                MPI_Send(&colMasterTemp_data[i][0], colElementCount, MPI_DOUBLE, i, 0, controlData.row_comm);
                std::cout << "Sent data" << std::endl;

                // Dont need to keep this data on the column master once its sent to the proper node!
                colMasterTemp_row[i].clear();
                colMasterTemp_col[i].clear();
                colMasterTemp_data[i].clear();
            }
        }

        // column masters recieve data from cluster master
        if (controlData.myId < controlData.clusterRows && controlData.myId > 0) {
            std::cout << controlData.myId << "column master reading data" << std::endl;
            MPI_Recv(&colElementCount, 1, MPI_INT, 0, 0, controlData.row_comm, MPI_STATUS_IGNORE);  // get number of elements
            MPI_Recv(&controlData.rowCount, 1, MPI_INT, 0, 0, controlData.row_comm, MPI_STATUS_IGNORE); // get number of rows (should be all)

            //std::cout << "rowCount = " << rowCount << ", colElementCount = " << colElementCount << std::endl;
            // resize the vectors to allow for the incoming data to be recieved
            std::vector <int> a,b;
            std::vector <double> c;
            colMasterTemp_row.push_back(a);
            colMasterTemp_col.push_back(b);
            colMasterTemp_data.push_back(c);
            //std::cout << "colMasterTemp_row[0].size() = " << colMasterTemp_row[0].size() << std::endl;
            //std::cout << "colMasterTemp_col[0].size() = " << colMasterTemp_col[0].size() << std::endl;
            //std::cout << "colMasterTemp_data[0].size() = " << colMasterTemp_data[0].size() << std::endl;

            colMasterTemp_row[0].resize(controlData.rowCount);
            colMasterTemp_col[0].resize(colElementCount);
            colMasterTemp_data[0].resize(colElementCount);

            //std::cout << "colMasterTemp_row[0].size() = " << colMasterTemp_row[0].size() << std::endl;
            //std::cout << "colMasterTemp_col[0].size() = " << colMasterTemp_col[0].size() << std::endl;
            //std::cout << "colMasterTemp_data[0].size() = " << colMasterTemp_data[0].size() << std::endl;
            //std::cout << "colMasterTemp_row.size() = " << colMasterTemp_row.size() << std::endl;

            MPI_Recv(&colMasterTemp_row[0][0], controlData.rowCount, MPI_INT, 0, 0, controlData.row_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&colMasterTemp_col[0][0], colElementCount, MPI_INT, 0, 0, controlData.row_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&colMasterTemp_data[0][0], colElementCount, MPI_DOUBLE, 0, 0, controlData.row_comm, MPI_STATUS_IGNORE);
        }
*/
        // Copy colMasterTemp data to local csr vectors for the column masters to work on after distrubting data to rows
        if (controlData.myId < controlData.clusterRows) {


        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Comm_free(&controlData.col_comm);
        MPI_Comm_free(&controlData.row_comm);
        MPI_Finalize();


    return 0;
}