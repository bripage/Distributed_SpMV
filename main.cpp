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
#include <unistd.h>


int main(int argc, char *argv[]) {
    controlData control;

    //Initialize the MPI environment
    MPI_Init(&argc, &argv);

    double overallStartTime = MPI_Wtime();

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &control.processCount);
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &control.myId);

    std::string argTemp;
    for (int i = 1; i < argc; i = i + 2) {
        argTemp = argv[i];
        if (argTemp == "-lm") {
            // load matrix from file.
            control.matrixFile = argv[i + 1];
        } else if (argTemp == "-lv") {
            // load dense vector from file.
            control.vectorFile = argv[i + 1];
        } else if (argTemp == "--master-only") {
            //run on MPI master only (no distribution, OpenMP as normal).
            std::string temp = argv[i + 1];
            if (temp == "true") {
                control.masterOnly = true;
            }
        } else if (argTemp == "--distribution-method") {
            //set number of OpenMP threads per node
            control.distributionMethod = argv[i + 1];
        } else if (argTemp == "--major-order") {
            std::string temp = argv[i + 1];
            if (temp == "col") {
                control.colMajor = true;
            }
        } else if (argTemp == "--omp-threads") {
            //set number of OpenMP threads per node
            control.ompThreads = atoi(argv[i + 1]);
        } else if (argTemp == "--help") {
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
    //
    // End MPI initialization
    //

    //std::vector<double> denseVector;

    csrSpMV masterData;
    if (!control.myId) {
        masterData.masterOnlySpMV(control);
    }


    //***********************************************//
    //     Select distribution method and actions    //
    //***********************************************//
    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<csrSpMV *> clusterColData;
    int colsPerColumn;

    double distributionStartTime = MPI_Wtime();

    if (!control.myId) {
        distribution_SplitMatrix(control, clusterColData);
    }

    double distributionEndTime = MPI_Wtime();

    // Create a pointer to the nodes csr object. Using pointers so that we do not have to copy any data on the
    // master node in order to assign it to this name
    csrSpMV *nodeCSR;

    if (control.myId == 0) {
        nodeCSR = clusterColData[0];
    } else {
        nodeCSR = new csrSpMV;
    }

    if (control.masterOnly != true) {
        // master to send data to cluster column masters
        if (control.myId == 0) {
            for (int i = 1; i < control.clusterCols; i++) {  // start at 1 since Master is the row master

	            //std::cout << "1 data size = " << clusterColData[0]->csrData.size() << std::endl;
                control.elementCount = clusterColData[i]->csrData.size();
	            //std::cout << "2 data size = " << clusterColData[0]->csrData.size() << std::endl;
                MPI_Send(&control.elementCount, 1, MPI_INT, i, 0, control.row_comm);
	            //std::cout << "3 data size = " << clusterColData[0]->csrData.size() << std::endl;
                MPI_Send(&control.rowCount, 1, MPI_INT, i, 0, control.row_comm);
	            //std::cout << "4 data size = " << clusterColData[0]->csrData.size() << std::endl;
                MPI_Send(&(clusterColData[i]->csrRows[0]), control.rowCount, MPI_INT, i, 0, control.row_comm);
	            //std::cout << "5 data size = " << clusterColData[0]->csrData.size() << std::endl;
                //MPI_Send(&(clusterColData[i]->csrCols[0]), clusterColData[i]->csrCols.size(), MPI_INT, i, 0,
                         control.row_comm);
	            //std::cout << "6 data size = " << clusterColData[0]->csrData.size() << std::endl;
                MPI_Send(&(clusterColData[i]->csrData[0]), clusterColData[i]->csrData.size(), MPI_DOUBLE, i, 0,
                         control.row_comm);
	            //std::cout << "7 data size = " << clusterColData[0]->csrData.size() << std::endl;
                MPI_Send(&(clusterColData[0]->denseVec[i * control.rowsPerNode]), control.rowsPerNode, MPI_DOUBLE, i, 0,
                         control.row_comm);
	            //std::cout << "8 data size = " << clusterColData[0]->csrData.size() << std::endl;
            }
	        //std::cout << "9 data size = " << clusterColData[0]->csrData.size() << std::endl;
            // delete and free column data that the master has already sent to the column masters, as it no longer needs
            // to be kept on the master
            for (int i = 1; i < clusterColData.size(); i++) {
                delete (clusterColData[i]);
            }
	        //std::cout << "10 data size = " << clusterColData[0]->csrData.size() << std::endl;
            clusterColData.erase(clusterColData.begin() + 1,
                                 clusterColData.end()); // remove pointers from clusterColData
	        //std::cout << "11 data size = " << clusterColData[0]->csrData.size() << std::endl;
        } else if (control.myId < control.clusterCols && control.myId != 0) {
            MPI_Recv(&control.elementCount, 1, MPI_INT, 0, 0, control.row_comm,
                     MPI_STATUS_IGNORE);  // get number of elements
            MPI_Recv(&control.rowCount, 1, MPI_INT, 0, 0, control.row_comm,
                     MPI_STATUS_IGNORE); // get number of rows (should be all)

            control.rowsPerNode = ceil(control.rowCount / (float) control.clusterRows);

            nodeCSR->csrRows.resize(control.rowCount);
            //nodeCSR->csrCols.resize(control.elementCount);
            nodeCSR->csrData.resize(control.elementCount);
            nodeCSR->denseVec.resize(control.rowsPerNode);

            MPI_Recv(&nodeCSR->csrRows[0], control.rowCount, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&nodeCSR->csrCols[0], control.elementCount, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&nodeCSR->csrData[0], control.elementCount, MPI_DOUBLE, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&nodeCSR->denseVec[0], control.rowsPerNode, MPI_DOUBLE, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
        }

        if (control.myId == 0) {
            nodeCSR->denseVec.erase(nodeCSR->denseVec.begin() + control.rowsPerNode, nodeCSR->denseVec.end());
        }
	    //std::cout << "12 data size = " << nodeCSR->csrData.size() << std::endl;
/*
	    usleep(100000*control.myId);
	    std::cout << std::endl;
	    std::cout << std::endl;
	    std::cout << "process " << control.myId << " rows: " << std::endl;
	    for (int j = 0; j < nodeCSR->csrRows.size(); j++) {
		    std::cout << nodeCSR->csrRows[j] << " ";
	    }
	    std::cout << std::endl;
	    std::cout << "process " << control.myId << " cols: " << std::endl;
	    for (int j = 0; j < nodeCSR->csrCols.size(); j++) {
		    std::cout << nodeCSR->csrCols[j] << " ";
	    }
	    std::cout << std::endl;
	    std::cout << "process " << control.myId << " data: " << std::endl;
	    for (int j = 0; j < nodeCSR->csrData.size(); j++) {
		    std::cout << nodeCSR->csrData[j] << " ";
	    }
	    std::cout << std::endl;
	    std::cout << "process " << control.myId << " denseVec: " << std::endl;
	    for (int j = 0; j < nodeCSR->denseVec.size(); j++) {
		    std::cout << nodeCSR->denseVec[j] << " ";
	    }
	    std::cout << std::endl;
	    std::cout << std::endl;
*/
        // column masters send data to row nodes
        if (control.myId < control.clusterCols) {
            for (int i = 1;
                 i < control.clusterRows; i++) {  // start at 1 since column master is the 0th node in the column
                int localRowCount;
                int firstElement = nodeCSR->csrRows[i * control.rowsPerNode];
                int lastElement;

                if (i == control.clusterRows - 1) {
                    lastElement = nodeCSR->csrData.size();
                    localRowCount = control.rowCount - (i * control.rowsPerNode);
                } else {
                    lastElement = nodeCSR->csrRows[(i + 1) * control.rowsPerNode];
                    localRowCount = control.rowsPerNode;
                }
                control.elementCount = (lastElement -
                                        firstElement); // how many elements the ith clusterCol has been assigned

                MPI_Send(&control.elementCount, 1, MPI_INT, i, 0, control.col_comm);
                MPI_Send(&localRowCount, 1, MPI_INT, i, 0, control.col_comm);
                MPI_Send(&control.rowsPerNode, 1, MPI_INT, i, 0, control.col_comm);

                MPI_Send(&(nodeCSR->csrRows[control.rowsPerNode * i]), localRowCount, MPI_INT, i, 0, control.col_comm);
                //MPI_Send(&(nodeCSR->csrCols[firstElement]), control.elementCount, MPI_INT, i, 0, control.col_comm);
                MPI_Send(&(nodeCSR->csrData[firstElement]), control.elementCount, MPI_DOUBLE, i, 0, control.col_comm);
                //MPI_Send(&(nodeCSR->denseVec[control.rowsPerNode * i]), localRowCount, MPI_DOUBLE, i, 0, control.col_comm);
            }

            // Erase the excess data on the column master that has already been distributed to its row nodes
            if (nodeCSR->csrRows.empty() || nodeCSR->csrCols.empty() || nodeCSR->csrData.empty()) {
                // empty, nothing to erase safely
            } else {
                int myLastData = nodeCSR->csrRows[control.rowsPerNode];
                nodeCSR->csrRows.erase(nodeCSR->csrRows.begin() + control.rowsPerNode, nodeCSR->csrRows.end());
                nodeCSR->csrCols.erase(nodeCSR->csrCols.begin() + myLastData, nodeCSR->csrCols.end());
                nodeCSR->csrData.erase(nodeCSR->csrData.begin() + myLastData, nodeCSR->csrData.end());
            }
        }

        if (control.myId < control.clusterCols) {
            nodeCSR->denseVec.erase(nodeCSR->denseVec.begin() + control.rowsPerNode, nodeCSR->denseVec.end());
        }

        if (control.myId >= control.clusterCols) {
            MPI_Recv(&control.elementCount, 1, MPI_INT, 0, 0, control.col_comm,
                     MPI_STATUS_IGNORE);  // get number of elements
            MPI_Recv(&control.rowCount, 1, MPI_INT, 0, 0, control.col_comm, MPI_STATUS_IGNORE); // get number of rows
            MPI_Recv(&control.rowsPerNode, 1, MPI_INT, 0, 0, control.col_comm, MPI_STATUS_IGNORE);

            nodeCSR->csrRows.resize(control.rowCount);
            //nodeCSR->csrCols.resize(control.elementCount);
            nodeCSR->csrData.resize(control.elementCount);
            nodeCSR->denseVec.resize(control.rowsPerNode);

            MPI_Recv(&nodeCSR->csrRows[0], control.rowCount, MPI_INT, 0, 0, control.col_comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&nodeCSR->csrCols[0], control.elementCount, MPI_INT, 0, 0, control.col_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&nodeCSR->csrData[0], control.elementCount, MPI_DOUBLE, 0, 0, control.col_comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&nodeCSR->denseVec[0], control.rowCount, MPI_DOUBLE, 0, 0, control.col_comm, MPI_STATUS_IGNORE);
        }

        // broadcast dense vector to column nodes
        MPI_Bcast(&nodeCSR->denseVec[0], control.rowsPerNode, MPI_DOUBLE, 0, control.col_comm);
    }

    std::vector<double> result;
    result.resize(control.rowsPerNode, 0.0);
/*
	usleep(100000*control.myId);
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "process " << control.myId << " rows: " << std::endl;
	for (int j = 0; j < nodeCSR->csrRows.size(); j++) {
		std::cout << nodeCSR->csrRows[j] << " ";
	}
	std::cout << std::endl;
	std::cout << "process " << control.myId << " cols: " << std::endl;
	for (int j = 0; j < nodeCSR->csrCols.size(); j++) {
		std::cout << nodeCSR->csrCols[j] << " ";
	}
	std::cout << std::endl;
	std::cout << "process " << control.myId << " data: " << std::endl;
	for (int j = 0; j < nodeCSR->csrData.size(); j++) {
		std::cout << nodeCSR->csrData[j] << " ";
	}
	std::cout << std::endl;
	std::cout << "process " << control.myId << " denseVec: " << std::endl;
	for (int j = 0; j < nodeCSR->denseVec.size(); j++) {
		std::cout << nodeCSR->denseVec[j] << " ";
	}
	std::cout << std::endl;
	std::cout << std::endl;
 */

    MPI_Barrier(MPI_COMM_WORLD);
    double spmvStartTime = MPI_Wtime();

    if (nodeCSR->csrData.size() > 0) {
        int ompThreadId, start, end, i, j, rowsPerThread;

        #pragma omp parallel num_threads(control.ompThreads) shared(nodeCSR, result) private(ompThreadId, start, end, i, j, rowsPerThread)
        {
            rowsPerThread = ceil(nodeCSR->csrRows.size() / control.ompThreads);
            ompThreadId = omp_get_thread_num();
            if (ompThreadId == control.ompThreads - 1) {
                for (i = ompThreadId * rowsPerThread; i < nodeCSR->csrRows.size(); i++) {
                    if (i == nodeCSR->csrRows.size() - 1) {
                        for (j = nodeCSR->csrRows[i] - nodeCSR->csrRows[0]; j < nodeCSR->csrData.size(); j++) {
                            result[i] += nodeCSR->csrData[j] * (double) nodeCSR->denseVec[i];
                        }
                    } else {
                        for (j = nodeCSR->csrRows[i] - nodeCSR->csrRows[0];
                             j < nodeCSR->csrRows[i + 1] - nodeCSR->csrRows[0]; j++) {
                            result[i] += nodeCSR->csrData[j] * (double) nodeCSR->denseVec[i];
                        }
                    }
                }
            } else {
                for (i = ompThreadId * rowsPerThread; i < (ompThreadId + 1) * rowsPerThread; i++) {
                    if (i == nodeCSR->csrRows.size() - 1) {
                        for (j = nodeCSR->csrRows[i] - nodeCSR->csrRows[0]; j < nodeCSR->csrData.size(); j++) {
                            result[i] += nodeCSR->csrData[j] * (double) nodeCSR->denseVec[i];
                        }
                    } else {
                        for (j = nodeCSR->csrRows[i] - nodeCSR->csrRows[0];
                             j < nodeCSR->csrRows[i + 1] - nodeCSR->csrRows[0]; j++) {
                            result[i] += nodeCSR->csrData[j] * (double) nodeCSR->denseVec[i];
                        }
                    }
                }
            }
        }
    }

    //MPI_Barrier(MPI_COMM_WORLD);
    double spmvEndTime = MPI_Wtime();
    double reductionStartTime, reductionEndTime;

    if (control.masterOnly != true) {
        if (control.myId == 0) {
            //nodeCSR->denseVec.clear();
            result.resize(control.rowsPerNode * control.clusterRows, 0.0);
        }
        reductionStartTime = MPI_Wtime();
        /*
         *      MPI REDUCE w/ SUM FROM COLUMNS TO ROW MASTER(S)
         */
        if (control.myId < control.clusterCols) {
            if (control.row_rank == 0) {
                MPI_Reduce(MPI_IN_PLACE, &result[0], control.rowsPerNode, MPI_DOUBLE, MPI_SUM, 0,
                           control.row_comm);
            } else {
                MPI_Reduce(&result[0], &result[0], control.rowsPerNode, MPI_DOUBLE, MPI_SUM, 0,
                           control.row_comm);
            }
        } else {
            if (control.row_rank == 0) {
                MPI_Reduce(MPI_IN_PLACE, &result[0], control.rowCount, MPI_DOUBLE, MPI_SUM, 0,
                           control.row_comm);
            } else {
                MPI_Reduce(&result[0], &result[0], control.rowCount, MPI_DOUBLE, MPI_SUM, 0,
                           control.row_comm);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        reductionEndTime = MPI_Wtime();

        /*
         *      MPI GATHER FROM ROW MASTERS TO GLOBAL MASTER
         */
        if (control.myId % control.clusterCols == 0) {
            if (control.col_rank == 0) {
                MPI_Gather(MPI_IN_PLACE, control.rowsPerNode, MPI_DOUBLE, &result[0], control.rowsPerNode,
                           MPI_DOUBLE, 0, control.col_comm);
            } else {
                MPI_Gather(&result[0], control.rowsPerNode, MPI_DOUBLE, &result[0], control.rowsPerNode,
                           MPI_DOUBLE, 0, control.col_comm);
            }
        }
    }


    MPI_Comm_free(&control.col_comm);
    MPI_Comm_free(&control.row_comm);
    double overallEndTime = MPI_Wtime();
    MPI_Finalize();


    //
    //                                      -- Testing / Verifcation --
    // This block is to be used in conjunction with single node, single thread/process SpMV and compares the results of
    // the "master only" SpMV against those of the distributed version's results.
    //
    if (control.myId == 0) {
        std::cout << std::endl;
        for (int i = 0; i < control.rowCount; i++) {
            //std::cout << "i = " << i << std::endl;
            if (std::abs(masterData.result[i] - result[i]) > 0.1) {
                //std::cout << "--- ERROR: result[" << i << "] DOES NOT MATCH ---" << std::endl;
	            std::cout << "row " << i << ": " << masterData.result[i] << " != " << result[i] << std::endl;
            }
        }
    }


    if (control.myId == 0) {
        //std::cout << std::endl << "Complete!" << std::endl;

        //double precision = MPI_Wtick();
        //std::cout << "precision = " << precision  << std::endl;
        //std::cout << "start = " << overallStartTime << ", end = " <<overallEndTime << std::endl;
        //std::cout << "Element Distribution Time: " << distributionEndTime - distributionStartTime << std::endl;
        //std::cout << "SmPV Time: " << spmvEndTime - spmvStartTime << std::endl;
        //std::cout << "Total time elapsed: " << overallEndTime - overallStartTime << std::endl;
        std::cout << distributionEndTime - distributionStartTime << "," << spmvEndTime - spmvStartTime << ","
                  << reductionEndTime - reductionStartTime << "," << overallEndTime - overallStartTime << std::endl;
    }


	return 0;
}
