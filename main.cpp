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
#include <unistd.h>
#include <sched.h>


int main(int argc, char *argv[]) {
    controlData control;
	double distributionEndTime, distributionStartTime, dataTransmissionEnd, dataTransmissionStart, spmvEndTime,
			spmvStartTime, reductionEndTime, reductionStartTime, masterGatherEnd, masterGatherStart, overallEndTime,
			overallStartTime;

	//Initialize the MPI environment
    MPI_Init(&argc, &argv);

    overallStartTime = MPI_Wtime();

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
	        std::string temp = argv[i + 1];
	        if (temp == "split") {
		        control.distributionMethod = 1;
	        } else if (temp == "balanced"){
		        control.distributionMethod = 2;
	        }
        } else if (argTemp == "--major-order") {
            std::string temp = argv[i + 1];
            if (temp == "col") {
                control.colMajor = true;
            }
        } else if (argTemp == "--omp-threads") {
            //set number of OpenMP threads per node
            control.ompThreads = atoi(argv[i + 1]);
        } else if (argTemp == "--use-barriers") {
	        //set number of OpenMP threads per node
	        std::string temp = argv[i + 1];
	        if (temp == "true") {
		        control.barrier = true;
	        }
	    } else if (argTemp == "-verify") {
	        //set number of OpenMP threads per node
	        std::string temp = argv[i + 1];
	        if (temp == "true") {
		        control.verify = true;
	        }
        } else if (argTemp == "-debug") {
	        //set number of OpenMP threads per node
	        std::string temp = argv[i + 1];
	        if (temp == "true") {
		        control.debug = true;
		        control.verify = true;
	        }
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

	// Get the name of the processor
	if (control.debug) {
		usleep(100000 * control.myId);
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		int name_len;
		MPI_Get_processor_name(processor_name, &name_len);
		std::cout << "Process " << control.myId << "is on process" << processor_name << std::endl;
	}

    //****************************************************//
    //  Create comm for each column/row of compute nodes  //
    //****************************************************//
	if (control.debug && control.myId == 0) std::cout << "Starting MPI Initialization" << std::endl;
	// get the number of MPI processes for work split and distribution
    control.clusterRows = sqrt(control.processCount);
    control.clusterCols = control.clusterRows; // works because cluster is expected to be square
    control.myCol = control.myId % control.clusterRows;
	//std::cout << "myId = " << control.myId << ", myCol = " << control.myCol << std::endl;

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
	if (control.debug && control.myId == 0) std::cout << "MPI Initialization complete" << std::endl;

    //***********************************************//
    //     Select distribution method and actions    //
    //***********************************************//
    std::vector<csrSpMV *> clusterColData;
    int colsPerColumn;

	if (control.debug && control.myId == 0) std::cout << "starting multinode distribution" << std::endl;
    distributionStartTime = MPI_Wtime();
    if (!control.myId) {
	    if (control.distributionMethod == 1) {
		    distribution_SplitMatrix(control, clusterColData);
	    } else if (control.distributionMethod == 2){
		    distribution_Balanced(control, clusterColData);
	    }
    }
	distributionEndTime = MPI_Wtime();
	if (control.debug && control.myId == 0) std::cout << "Multinode distribution complete" << std::endl;

    //
    // Find Distribution of Non-Zeros on both sequential and distributed distributions
    //
	if (control.debug && control.myId == 0) std::cout << "Starting verification" << std::endl;
	csrSpMV masterData;
	if (control.verify) {
        if (!control.myId) {
            std::vector<int> seqDist(control.rowCount, 0); // sequential distribution
            std::vector<int> distDist(control.rowCount, 0); //distributed distribution

            masterData.masterOnlySpMV(control, seqDist); // perform sequential SpMV on master process only
	        if (control.debug && control.myId == 0) std::cout << "Sequential SpMV complete" << std::endl;

            // Get Distributed Row NNZ Count
            for (int i = 0; i < control.clusterCols; i++){
                for (int j = 0; j < control.rowCount; j++){
                    if(j == control.rowCount-1) {
                        distDist[j] += clusterColData[i]->csrData.size() - clusterColData[i]->csrRows[j];
                    } else {
                        distDist[j] += clusterColData[i]->csrRows[j+1] - clusterColData[i]->csrRows[j];
                    }
                }
            }

            int totalInvalidRows = 0;
            for (int i = 0; i < control.rowCount; i++){
                if (seqDist[i] != distDist[i]){
                    totalInvalidRows++;
                }
            }
            if (totalInvalidRows != 0){
                std::cout << "Total Incorrectly Distributed Rows = " << totalInvalidRows << std::endl;
            }
        }
    }
	//
    // End distribution calculation
	if (control.debug && control.myId == 0) std::cout << "Verification complete" << std::endl;


    // Create a pointer to the nodes csr object. Using pointers so that we do not have to copy any data on the
    // master node in order to assign it to this name
    csrSpMV *nodeCSR;
	std::vector <double> result;

    if (control.myId == 0) {
        nodeCSR = clusterColData[0];
    } else {
        nodeCSR = new csrSpMV;
    }

	if (control.barrier && control.myRow == 0) MPI_Barrier(control.row_comm);
	dataTransmissionStart = MPI_Wtime();
	if (control.distributionMethod == 1) {
		if (control.masterOnly != true) {
			if (control.debug && control.myId == 0) std::cout << "Sending data to column masters" << std::endl;
			// master to send data to cluster column masters
			if (control.myId == 0) {
				for (int i = 1; i < control.clusterCols; i++) {  // start at 1 since Master is the row master
					control.elementCount = clusterColData[i]->csrData.size();
					MPI_Send(&control.elementCount, 1, MPI_INT, i, 0, control.row_comm);
					MPI_Send(&control.rowCount, 1, MPI_INT, i, 0, control.row_comm);
					MPI_Send(&(clusterColData[i]->csrRows[0]), control.rowCount, MPI_INT, i, 0, control.row_comm);
					MPI_Send(&(clusterColData[i]->csrCols[0]), clusterColData[i]->csrCols.size(), MPI_INT, i, 0,
					         control.row_comm);
					MPI_Send(&(clusterColData[i]->csrData[0]), clusterColData[i]->csrData.size(), MPI_DOUBLE, i, 0,
					         control.row_comm);
					MPI_Send(&(clusterColData[0]->denseVec[i * control.rowsPerNode]), control.rowsPerNode, MPI_DOUBLE,
					         i, 0,
					         control.row_comm);
				}
				// delete and free column data that the master has already sent to the column masters, as it no longer needs
				// to be kept on the master
				for (int i = 1; i < clusterColData.size(); i++) {
					delete (clusterColData[i]);
				}
				clusterColData.erase(clusterColData.begin() + 1, clusterColData.end()); //remove transmitted data
				nodeCSR->denseVec.erase(nodeCSR->denseVec.begin() + control.rowsPerNode, nodeCSR->denseVec.end());
			} else if (control.myId < control.clusterCols && control.myId != 0) {
				MPI_Recv(&control.elementCount, 1, MPI_INT, 0, 0, control.row_comm,
				         MPI_STATUS_IGNORE);  // get number of elements
				MPI_Recv(&control.rowCount, 1, MPI_INT, 0, 0, control.row_comm,
				         MPI_STATUS_IGNORE); // get number of rows (should be all)

				control.rowsPerNode = ceil(control.rowCount / (float) control.clusterRows);
				control.colsPerNode = control.rowsPerNode;

				nodeCSR->csrRows.resize(control.rowCount);
				nodeCSR->csrCols.resize(control.elementCount);
				nodeCSR->csrData.resize(control.elementCount);
				nodeCSR->denseVec.resize(control.rowsPerNode);

				MPI_Recv(&nodeCSR->csrRows[0], control.rowCount, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrCols[0], control.elementCount, MPI_INT, 0, 0, control.row_comm,
				         MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrData[0], control.elementCount, MPI_DOUBLE, 0, 0, control.row_comm,
				         MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->denseVec[0], control.rowsPerNode, MPI_DOUBLE, 0, 0, control.row_comm,
				         MPI_STATUS_IGNORE);
			}
			if (control.debug && control.myId == 0) std::cout << "Sending to column masters complete" << std::endl;

			// column masters send data to row nodes
			if (control.barrier) MPI_Barrier(control.col_comm);
			if (control.myId < control.clusterCols) {
				if (control.debug && control.myId == 0)
					std::cout << "Column Masters Sending Data to col members" << std::endl;
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

					MPI_Send(&(nodeCSR->csrRows[control.rowsPerNode * i]), localRowCount, MPI_INT, i, 0,
					         control.col_comm);
					MPI_Send(&(nodeCSR->csrCols[firstElement]), control.elementCount, MPI_INT, i, 0, control.col_comm);
					MPI_Send(&(nodeCSR->csrData[firstElement]), control.elementCount, MPI_DOUBLE, i, 0,
					         control.col_comm);
					//MPI_Send(&(nodeCSR->denseVec[control.rowsPerNode * i]), localRowCount, MPI_DOUBLE, i, 0, control.col_comm);
				}

				// Erase the excess data on the column master that has already been distributed to its row nodes
				int myLastData = nodeCSR->csrRows[control.rowsPerNode];
				if (!(nodeCSR->csrRows.empty())) {
					nodeCSR->csrRows.erase(nodeCSR->csrRows.begin() + control.rowsPerNode, nodeCSR->csrRows.end());
				}
				if (!(nodeCSR->csrCols.empty())) {
					nodeCSR->csrCols.erase(nodeCSR->csrCols.begin() + myLastData, nodeCSR->csrCols.end());
				}
				if (!(nodeCSR->csrData.empty())) {
					nodeCSR->csrData.erase(nodeCSR->csrData.begin() + myLastData, nodeCSR->csrData.end());
				}
				nodeCSR->rebase(control.myCol * control.colsPerNode);

				nodeCSR->denseVec.erase(nodeCSR->denseVec.begin() + control.rowsPerNode, nodeCSR->denseVec.end());
			} else if (control.myId >= control.clusterCols) {
				MPI_Recv(&control.elementCount, 1, MPI_INT, 0, 0, control.col_comm,
				         MPI_STATUS_IGNORE);  // get number of elements
				MPI_Recv(&control.rowCount, 1, MPI_INT, 0, 0, control.col_comm,
				         MPI_STATUS_IGNORE); // get number of rows
				MPI_Recv(&control.rowsPerNode, 1, MPI_INT, 0, 0, control.col_comm, MPI_STATUS_IGNORE);

				nodeCSR->csrRows.resize(control.rowCount);
				nodeCSR->csrCols.resize(control.elementCount);
				nodeCSR->csrData.resize(control.elementCount);
				nodeCSR->denseVec.resize(control.rowsPerNode);

				MPI_Recv(&nodeCSR->csrRows[0], control.rowCount, MPI_INT, 0, 0, control.col_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrCols[0], control.elementCount, MPI_INT, 0, 0, control.col_comm,
				         MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrData[0], control.elementCount, MPI_DOUBLE, 0, 0, control.col_comm,
				         MPI_STATUS_IGNORE);
				//MPI_Recv(&nodeCSR->denseVec[0], control.rowCount, MPI_DOUBLE, 0, 0, control.col_comm, MPI_STATUS_IGNORE);

				control.colsPerNode = control.rowsPerNode;
				nodeCSR->rebase(control.myCol * control.colsPerNode);
			}

			// broadcast dense vector to column nodes
			MPI_Bcast(&nodeCSR->denseVec[0], control.rowsPerNode, MPI_DOUBLE, 0, control.col_comm);
		}
		dataTransmissionEnd = MPI_Wtime();
		if (control.debug && control.myId == 0) {
			std::cout << "Col master sending to col members complete" << std::endl;
			usleep(10000000 * control.myId);
		}
		
		result.resize(control.rowsPerNode, 0.0);

		if (control.barrier) MPI_Barrier(MPI_COMM_WORLD);
		spmvStartTime = MPI_Wtime();

		if (control.debug && control.myId == 0) std::cout << "Starting SpMV computation" << std::endl;
		if (nodeCSR->csrData.size() > 0) {
			int ompThreadId, ompCPUId, start, end, i, j, k, rowsPerThread, rowEnd;

#pragma omp parallel num_threads(control.ompThreads) shared(nodeCSR, result) private(ompThreadId, ompCPUId, start, end, i, j, k, rowsPerThread, rowEnd)
			{
				rowsPerThread = ceil(nodeCSR->csrRows.size() / (double) control.ompThreads);
				ompThreadId = omp_get_thread_num();
				ompCPUId = sched_getcpu();

				if (control.debug) {
					usleep(100000 * ompThreadId);
					std::cout << "Rank " << control.myId << ", Thread " << ompThreadId << " on core " << ompCPUId
					          << std::endl;
				}


				if (ompThreadId == control.ompThreads - 1) {
					rowEnd = nodeCSR->csrRows.size();
				} else {
					rowEnd = (ompThreadId + 1) * rowsPerThread;
				}

				if (ompThreadId == control.ompThreads - 1) {
					for (i = ompThreadId * rowsPerThread; i < nodeCSR->csrRows.size(); i++) {
						if (i == nodeCSR->csrRows.size() - 1) {
							for (j = nodeCSR->csrRows[i]; j < nodeCSR->csrData.size(); j++) {
								result[i] += nodeCSR->csrData[j] * (double) nodeCSR->denseVec[nodeCSR->csrCols[j]];
							}
						} else {
							for (j = nodeCSR->csrRows[i]; j < nodeCSR->csrRows[i + 1]; j++) {
								result[i] += nodeCSR->csrData[j] * (double) nodeCSR->denseVec[nodeCSR->csrCols[j]];
							}
						}
					}
				} else {
					for (k = ompThreadId * rowsPerThread; k < rowEnd; k++) {
						for (j = nodeCSR->csrRows[k]; j < nodeCSR->csrRows[k + 1]; j++) {
							result[k] += nodeCSR->csrData[j] * (double) nodeCSR->denseVec[nodeCSR->csrCols[j]];
						}
					}
				}
			}
		}
		if (control.barrier) MPI_Barrier(MPI_COMM_WORLD);
		if (control.debug && control.myId == 0) std::cout << "SpMV computation complete" << std::endl;
		double spmvEndTime = MPI_Wtime();

		if (control.masterOnly != true) {
			if (control.myId == 0) {
				result.resize(control.rowsPerNode * control.clusterRows, 0.0);
			}
			reductionStartTime = MPI_Wtime();
			/*
			 *      MPI REDUCE w/ SUM FROM COLUMNS TO ROW MASTER(S)
			 */
			if (control.debug && control.myId == 0) std::cout << "Starting MPI Reduction" << std::endl;
			if (control.myId < control.clusterCols) {
				if (control.row_rank == 0) {
					MPI_Reduce(MPI_IN_PLACE, &result[0], control.rowsPerNode, MPI_DOUBLE, MPI_SUM, 0, control.row_comm);
				} else {
					MPI_Reduce(&result[0], &result[0], control.rowsPerNode, MPI_DOUBLE, MPI_SUM, 0, control.row_comm);
				}
			} else {
				if (control.row_rank == 0) {
					MPI_Reduce(MPI_IN_PLACE, &result[0], control.rowCount, MPI_DOUBLE, MPI_SUM, 0, control.row_comm);
				} else {
					MPI_Reduce(&result[0], &result[0], control.rowCount, MPI_DOUBLE, MPI_SUM, 0, control.row_comm);
				}
			}
			if (control.barrier) MPI_Barrier(MPI_COMM_WORLD);
			if (control.debug && control.myId == 0) std::cout << "MPI Reduction complete" << std::endl;
			reductionEndTime = MPI_Wtime();

			/*
			 *      MPI GATHER FROM ROW MASTERS TO GLOBAL MASTER
			 */
			if (control.barrier) MPI_Barrier(control.col_comm);
			masterGatherStart = MPI_Wtime();
			if (control.debug && control.myId == 0) std::cout << "Starting MPI Gather" << std::endl;
			if (control.myId % control.clusterCols == 0) {
				if (control.col_rank == 0) {
					MPI_Gather(MPI_IN_PLACE, control.rowsPerNode, MPI_DOUBLE, &result[0], control.rowsPerNode,
					           MPI_DOUBLE, 0, control.col_comm);
				} else {
					MPI_Gather(&result[0], control.rowsPerNode, MPI_DOUBLE, &result[0], control.rowsPerNode,
					           MPI_DOUBLE, 0, control.col_comm);
				}
			}
			if (control.debug && control.myId == 0) std::cout << "MPI Gather complete" << std::endl;
			masterGatherEnd = MPI_Wtime();
		}
	}

	if (control.debug && control.myId == 0) std::cout << "Starting Finalization" << std::endl;
    MPI_Comm_free(&control.col_comm);
    MPI_Comm_free(&control.row_comm);
    overallEndTime = MPI_Wtime();
    MPI_Finalize();
	if (control.debug && control.myId == 0) std::cout << "Finalization complete" << std::endl;

    //
    //                                      -- Testing / Verifcation --
    // This block is to be used in conjunction with single node, single thread/process SpMV and compares the results of
    // the "master only" SpMV against those of the distributed version's results.
    //
	if (control.debug && control.myId == 0) std::cout << "Starting Verification" << std::endl;
	if (control.verify) {
        if (control.myId == 0) {
            std::cout << std::endl;
            int incorrectRowCount = 0;
            for (int i = 0; i < control.rowCount; i++) {
                //std::cout << "i = " << i << std::endl;
                if (std::abs(masterData.result[i] - result[i]) > 0.1) {
                    //std::cout << "--- ERROR: result[" << i << "] DOES NOT MATCH ---" << std::endl;
                    if (incorrectRowCount < 50) {
                        std::cout << "row " << i << ": " << masterData.result[i] << " != " << result[i] << std::endl;
                    }
                    incorrectRowCount++;
                }
            }
            std::cout << "Incorecct Rows: " << incorrectRowCount << std::endl;
        }
    }
	if (control.debug && control.myId == 0) std::cout << "Verification complete" << std::endl;

	if (control.debug && control.myId == 0) std::cout << "Starting Output" << std::endl;
    if (control.myId == 0) {
	    if (control.barrier == false) {
		    std::cout << distributionEndTime - distributionStartTime << "," << overallEndTime - overallStartTime << std::endl;
	    } else {
		    std::cout << distributionEndTime - distributionStartTime << ","
		              << dataTransmissionEnd - dataTransmissionStart << "," << spmvEndTime - spmvStartTime << ","
		              << reductionEndTime - reductionStartTime << "," << masterGatherEnd - masterGatherStart << ","
		              << overallEndTime - overallStartTime << std::endl;
	    }

    }
	if (control.debug && control.myId == 0) std::cout << "Output complete" << std::endl;


	return 0;
}
