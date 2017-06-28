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
    controlData control;
    std::string argTemp;

    for (int i = 1; i < argc; i= i+2){
        argTemp = argv[i];
        if (argTemp == "-lm"){
            // load matrix from file.
            control.matrixFile = argv[i+1];
        } else if (argTemp == "-lv"){
            // load dense vector from file.
            control.vectorFile = argv[i+1];
        } else if (argTemp == "-s"){
            //run on MPI master only (no distribution, OpenMP as normal).
            if (argv[i+1] == "true"){
                control.masterOnly = true;
            }
        } else if (argTemp == "--distribution-method"){
            //set number of OpenMP threads per node
            control.distributionMethod = argv[i+1];
        } else if (argTemp == "--major-order") {
            std::string temp = argv[i + 1];
            if (temp == "col") {
                control.colMajor = true;
            }
        } else if (argTemp == "--omp-threads"){
                //set number of OpenMP threads per node
            control.ompThreads = atoi(argv[i+1]);
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
    MPI_Comm_size(MPI_COMM_WORLD, &control.processCount);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &control.myId);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Verify MPI initialised properly
    if (control.processCount == 1 && control.masterOnly == false) {
        printf("MPI Initialization failed -- KILLING!");
        MPI_Finalize();
        exit(0);
    }

    //****************************************************//
    //  Create comm for each column/row of compute nodes  //
    //****************************************************//
    // get the number of MPI processes for work split and distribution
    control.clusterRows = sqrt(control.processCount);
    control.clusterCols = control.clusterRows; // works because cluster is expected to be square
    control.myCol = control.myId % control.clusterRows;

    // Create a duplicate of MPI_COMM_WORLD that can be used to split
    MPI_Comm dupCommWorldCol;
    MPI_Comm_dup(MPI_COMM_WORLD, &dupCommWorldCol);

    // Split dupCommWorld into comm's for each node column, based on each processes original myId
    //control.row_comm = new MPI_Comm;
    MPI_Comm_split(dupCommWorldCol, control.myCol, control.myId, &control.col_comm);

    //int col_rank, col_size;
    MPI_Comm_rank(control.col_comm, &control.col_rank);
    MPI_Comm_size(control.col_comm, &control.col_size);

    // get the number of MPI processes for work split and distribution
    control.myRow = control.myId / control.clusterRows;
    // Create a duplicate of MPI_COMM_WORLD that can be used to split
    MPI_Comm dupCommWorldRow;
    MPI_Comm_dup(MPI_COMM_WORLD, &dupCommWorldRow);

    // Split dupCommWorld into comm's for each node column, based on each processes original myId
    //control.row_comm = new MPI_Comm;
    MPI_Comm_split(dupCommWorldRow, control.myRow, control.myId, &control.row_comm);

    MPI_Comm_rank(control.row_comm, &control.row_rank);
    MPI_Comm_size(control.row_comm, &control.row_size);

    // create vectors to hold sparse matrix data once converted to CSR format
    std::vector<double> denseVector;


    csrSpMV masterData;
    if (!control.myId) {
        masterData.masterOnlySpMV(control);
    }


    //***********************************************//
    //     Select distribution method and actions    //
    //***********************************************//
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<csrSpMV*> clusterColData;
    int colsPerColumn, colElementCount;

    if (!control.myId) {
        distribution_SplitMatrix(control, clusterColData);

    }

    // create a pointer to the nodes csr object
    // using pointers so that we do not have to copy any data on the master node in order to assign it to this name
    csrSpMV* nodeCSR;

    if (!control.myId){
        nodeCSR = clusterColData[0];
    } else {
        nodeCSR = new csrSpMV;
    }

    // master to send data to cluster column masters
    if (!control.myId){
        for (int i = 1; i < control.clusterCols; i++){  // start at 1 since Master is the row master
            /*
            int firstElement = clusterColData[i]->csrRows[0];
            int lastElement;
            if (i == control.clusterCols -1){
                lastElement = clusterColData.size();
            } else {
                lastElement = clusterColData[i]->csrRows[(i+1) * control.rowsPerNode] - 1;
            }

            colElementCount = lastElement - firstElement; // how many elements the ith clusterCol has been assigned
            */

            colElementCount = clusterColData[i]->csrData.size();
            MPI_Send(&colElementCount, 1, MPI_INT, i, 0, control.row_comm);
            std::cout << "Sent colElementCount" << std::endl;
            MPI_Send(&control.rowCount, 1, MPI_INT, i, 0, control.row_comm);
            std::cout << "Sent rowCount" << std::endl;
            MPI_Send(&(clusterColData[i]->csrRows[0]), control.rowCount, MPI_INT, i, 0, control.row_comm);
            std::cout << "Sent rows" << std::endl;
            MPI_Send(&(clusterColData[i]->csrCols[0]), clusterColData[i]->csrCols.size(), MPI_INT, i, 0, control.row_comm);
            std::cout << "Sent cols" << std::endl;
            MPI_Send(&(clusterColData[i]->csrData[0]), clusterColData[i]->csrData.size(), MPI_DOUBLE, i, 0, control.row_comm);
            std::cout << "Sent data" << std::endl;
        }

        // delete and free column data that the master has already sent to the column masters, as it no longer needs
        // to be kept on the master
        for (int i = 1; i < clusterColData.size(); i++){
            delete(clusterColData[i]);
        }
        clusterColData.erase(clusterColData.begin()+1, clusterColData.end()); // remove pointers from clusterColData
    } else if (control.myId < control.clusterCols){
        std::cout << control.myId << "column master reading data" << std::endl;
        MPI_Recv(&colElementCount, 1, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);  // get number of elements
        MPI_Recv(&control.rowCount, 1, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE); // get number of rows (should be all)

        // resize the vectors to allow for the incoming data to be recieved
        nodeCSR->csrRows.resize(control.rowCount);
        nodeCSR->csrCols.resize(colElementCount);
        nodeCSR->csrData.resize(colElementCount);

        MPI_Recv(&nodeCSR->csrRows[0], control.rowCount, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
        MPI_Recv(&nodeCSR->csrCols[0], colElementCount, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
        MPI_Recv(&nodeCSR->csrData[0], colElementCount, MPI_DOUBLE, 0, 0, control.row_comm, MPI_STATUS_IGNORE);

        control.rowsPerNode = ceil(control.rowCount / (float)control.clusterRows);
    }

    // column masters send data to row nodes
    if (control.myId < control.clusterCols){
        for (int i = 1; i < control.clusterCols; i++){  // start at 1 since column master is the 0th node in the column

        }



        // erase the excess data on the column master that has already been distributed to its row nodes
        if (nodeCSR->csrRows.empty() || nodeCSR->csrCols.empty() || nodeCSR->csrData.empty()){
            std::cout << " empty! " << std::endl;
        } else {
            std::cout << "control.rowsPerNode = " << control.rowsPerNode << std::endl;
            std::cout << " nodeCSR->csrRows[control.rowsPerNode] = " << nodeCSR->csrRows[control.rowsPerNode] << std::endl;
            int myLastData = nodeCSR->csrRows[control.rowsPerNode] - 1;
            std::cout << "myLastData = " << myLastData << std::endl;
            nodeCSR->csrRows.erase(nodeCSR->csrRows.begin() + control.rowsPerNode, nodeCSR->csrRows.end());
            nodeCSR->csrCols.erase(nodeCSR->csrCols.begin() + myLastData, nodeCSR->csrCols.end());
            nodeCSR->csrData.erase(nodeCSR->csrData.begin() + myLastData, nodeCSR->csrData.end());
        }
    }




    //nodeCSR->nodeSpMV(control);

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Comm_free(&control.col_comm);
        MPI_Comm_free(&control.row_comm);
        MPI_Finalize();


    return 0;
}