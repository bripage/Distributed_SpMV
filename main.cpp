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
	        } else if (temp == "bv2"){
		        control.distributionMethod = 3;
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
        } else if (argTemp == "--show-matrix-info") {
	        //set number of OpenMP threads per node
	        std::string temp = argv[i + 1];
	        if (temp == "true") {
		        control.matrixInfo = true;
	        }
        }else if (argTemp == "--show-distribution") {
		    //set number of OpenMP threads per node
		    std::string temp = argv[i + 1];
		    if (temp == "true") {
		        control.distributionInfo = true;
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
    if (control.processCount == 1 && control.masterOnly == false) control.masterOnly = true; //you meant to put true ;)

	// Get the name of the processor
	if (control.debug) {
		//usleep(100000 * control.myId);
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
	if (control.debug && control.myId == 0) std::cout << "clusterCols = " << control.clusterCols << ", clusterRows = " << control.clusterRows << std::endl;
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
	std::vector <double> gatheredResult;
	std::vector <int> displacements, rowCounts;
    int colsPerColumn;

	if (control.debug && control.myId == 0) std::cout << "starting multinode distribution" << std::endl;
    distributionStartTime = MPI_Wtime();
    if (!control.myId) {
	    if (control.distributionMethod == 1) {
		    distribution_SplitMatrix(control, clusterColData);
	    } else if (control.distributionMethod == 2){
		    distribution_Balanced(control, clusterColData);
	    } else if (control.distributionMethod == 3) {
		    //std::cout << "inside distribution" << std::endl;
		    distribution_Balanced(control, clusterColData);

		    //displacements.resize(control.processCount,0);
		    //rowCounts.resize(control.processCount,0);
		    int rowsToGather = 0;

		    //std::cout << "assignedRows[0] = " << clusterColData[0]->assignedRowIds.size() << std::endl;
		    //std::cout << "assignedRows[1] = " << clusterColData[1]->assignedRowIds.size() << std::endl;

		    //rowCounts.push_back(clusterColData[0]->processData[1]);
		    for (int i = 0; i < control.processCount; i++){
			    //std::cout << "i = " << i << std::endl;
			    rowCounts.push_back(clusterColData[i%control.clusterCols]->processData[((i/control.clusterRows)*2)+1]);
			    //std::cout << "i = " << i << ", rowcounts[" << i << "] = " << rowCounts[i] << std::endl;
			    rowsToGather += clusterColData[i%control.clusterCols]->processData[((i/control.clusterRows)*2)+1];
			    if (i == 0){
				    displacements.push_back(0);
			    } else {
				    displacements.push_back(displacements[displacements.size() - 1] + rowCounts[i - 1]);
			    }

			    for (int j = 0; j < rowCounts[i]; j++){
				    control.rowDistribution.push_back(clusterColData[i%control.clusterCols]->assignedRowIds[clusterColData[i%control.clusterCols]->processData[((i/control.clusterRows)*2)+1]+j]);
			        //std::cout << "control.rowDistribution[" << i << "] = "
			        //          << clusterColData[i%control.clusterCols]->assignedRowIds[clusterColData[i%control.clusterCols]->processData[((i/control.clusterRows)*2)+1]+j]
			        //          << std::endl;
			    }
		    }

		    for (int i = 0; i < displacements.size(); i++){
			    //std::cout << "displacement[" << i << "] = " << displacements[i] << ", rowcounts[" << i << "] = " << rowCounts[i] << std::endl;
		    }
		    
		    //gatheredResult.resize(rowsToGather, 0.0);
	    }
    }
	distributionEndTime = MPI_Wtime();
	if (control.debug && control.myId == 0) std::cout << "Multinode distribution complete" << std::endl;

    //
    // Find Distribution of Non-Zeros on both sequential and distributed distributions
	// This checks to make sure that all nnz have been assigned and that none are forgotten or used twice
	//
	if (control.debug && control.myId == 0) std::cout << "Starting verification" << std::endl;
	csrSpMV masterData;
	if (control.verify) {
        if (!control.myId) {
	        if (control.distributionMethod == 1) {
		        std::vector<int> seqDist(control.rowCount, 0); // sequential distribution
		        std::vector<int> distDist(control.rowCount, 0); //distributed distribution

		        masterData.masterOnlySpMV(control, seqDist); // perform sequential SpMV on master process only
		        if (control.debug && control.myId == 0) std::cout << "Sequential SpMV complete" << std::endl;

		        // Get Distributed Row NNZ Count
		        for (int i = 0; i < control.clusterCols; i++) {
			        for (int j = 0; j < control.rowCount; j++) {
				        if (j == control.rowCount - 1) {
					        distDist[j] += clusterColData[i]->csrData.size() - clusterColData[i]->csrRows[j];
				        } else {
					        distDist[j] += clusterColData[i]->csrRows[j + 1] - clusterColData[i]->csrRows[j];
				        }
			        }
		        }

		        int totalInvalidRows = 0;
		        for (int i = 0; i < control.rowCount; i++) {
			        if (seqDist[i] != distDist[i]) {
				        totalInvalidRows++;
			        }
		        }
		        if (totalInvalidRows != 0) {
			        std::cout << "Total Incorrectly Distributed Rows = " << totalInvalidRows << std::endl;
		        }
	        } else {
		        // verification for balance distribution
		        std::vector<int> seqDist(control.rowCount, 0); // sequential distribution
		        std::cout << "made the seqDist vector" << std::endl;
		        masterData.masterOnlySpMV(control, seqDist); // perform sequential SpMV on master process only
	        }
        }
    }
	//
    // End distribution calculation
	if (control.debug && control.myId == 0) std::cout << "Verification complete" << std::endl;

    if (control.matrixInfo && control.myId == 0) {
	    if (control.distributionMethod == 1) {
		    int seqProcSum = 0;
		    double seqProcAvg, seqProcAvgDiff = 0.0, seqStandardDeviation;
		    std::vector<int> seqDist(control.rowCount, 0); // sequential distribution

		    for (int i = 0; i < control.processCount; i++) {
			    int firstElement = clusterColData[i%control.clusterCols]->csrRows[(i/control.clusterRows) * control.rowsPerNode];
			    int lastElement, nnzCount;

			    if (i/control.clusterRows == control.clusterRows - 1) {
				    lastElement = clusterColData[i%control.clusterCols]->csrData.size();
			    } else {
				    lastElement = clusterColData[i%control.clusterCols]->csrRows[((i/control.clusterRows) + 1) * control.rowsPerNode];
			    }
			    nnzCount = (lastElement - firstElement);
			    seqDist[i] = nnzCount;
			    seqProcSum += nnzCount;
			    std::cout << "Process " << i << ": " << seqDist[i] << std::endl;
		    }

		    seqProcAvg = seqProcSum / (double) control.processCount;
		    for (int i = 0; i < control.processCount; i++) {
			    seqProcAvgDiff += (seqDist[i] - seqProcAvg) * (seqDist[i] - seqProcAvg);
		    }
		    seqStandardDeviation = sqrt((1.0 / control.processCount) * seqProcAvgDiff);
		    std::cout << "NNZ Per Process Standard Deviation = " << seqStandardDeviation << std::endl;

		    if (control.debug && control.myId == 0)
			    std::cout << "Done Determining NNZ Per Process Standard Deviation" << std::endl;
	    } else if (control.distributionMethod == 2) {
		    //get distribution process averages
		    if (control.debug && control.myId == 0)
			    std::cout << "Starting Determining NNZ Per Process Standard Deviation" << std::endl;

		    int distProcSum = 0;
		    double distProcAvg, distProcAvgDiff = 0.0, distStandardDeviation;
		    std::vector<int> distDist(control.processCount, 0); //distributed distribution

		    for (int i = 0; i < control.processCount; i++) {
			    distDist[i] = clusterColData[i % control.clusterCols]->processData[(i / control.clusterRows) * 2];
			    distProcSum += distDist[i];
			    std::cout << "Process " << i << ": " << distDist[i] << std::endl;
		    }
		    distProcAvg = distProcSum / (double) control.processCount;
		    for (int i = 0; i < control.processCount; i++) {
			    distProcAvgDiff += (distDist[i] - distProcAvg) * (distDist[i] - distProcAvg);
		    }
		    distStandardDeviation = sqrt((1.0 / control.processCount) * distProcAvgDiff);
		    std::cout << "NNZ Per Process Standard Deviation = " << distStandardDeviation << std::endl;

		    if (control.debug && control.myId == 0)
			    std::cout << "Done Determining NNZ Per Process Standard Deviation" << std::endl;
	    }
    }

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
			//usleep(10000000 * control.myId);
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
					//usleep(100000 * ompThreadId);
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
		spmvEndTime = MPI_Wtime();

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
	} else if (control.distributionMethod == 2){    // balanced distribution
		if (control.masterOnly != true) {
			if (control.debug && control.myId == 0) std::cout << "Sending data to column masters" << std::endl;
			// master to send data to cluster column masters
			if (control.myId == 0) {
				for (int i = 1; i < control.clusterCols; i++) {  // start at 1 since Master is the row master
					MPI_Send(&control.rowCount, 1, MPI_INT, i, 0, control.row_comm);
					//std::cout << "preparing to send " << clusterColData[i]->processData.size() << " processData elements to " << i << ": ";
					//for (int j = 0; j < clusterColData[i]->processData.size(); j++){
					//	std::cout << "for  = " << j << ": ";
					//	std::cout << clusterColData[i]->processData[j] << std::endl;
					//}
					//std::cout << std::endl;
					//std::cout  << "sending processData to " << i << std::endl;
					MPI_Send(&(clusterColData[i]->processData[0]), control.clusterRows*2, MPI_INT, i, 0, control.row_comm);
					MPI_Send(&(clusterColData[i]->assignedRowIds[0]), clusterColData[i]->csrRows.size(), MPI_INT, i, 0, control.row_comm);
					//std::cout << "sending csrRows to " << i << std::endl;
					MPI_Send(&(clusterColData[i]->csrRows[0]), clusterColData[i]->csrRows.size(), MPI_INT, i, 0,
					         control.row_comm);
					//std::cout << "sending " << clusterColData[i]->csrCols.size() << " csrCols to " << i << std::endl;
					MPI_Send(&(clusterColData[i]->csrCols[0]), clusterColData[i]->csrCols.size(), MPI_INT, i, 0,
					         control.row_comm);
					//std::cout << "sending " << clusterColData[i]->csrData.size() << " csrData to " << i << std::endl;
					MPI_Send(&(clusterColData[i]->csrData[0]), clusterColData[i]->csrData.size(), MPI_DOUBLE, i, 0,
					         control.row_comm);
				}
				// delete and free column data that the master has already sent to the column masters, as it no longer needs
				// to be kept on the master
				for (int i = 1; i < clusterColData.size(); i++) {
					delete (clusterColData[i]);
				}
				//clusterColData.erase(clusterColData.begin() + 1, clusterColData.end()); //remove transmitted data
				//nodeCSR->denseVec.erase(nodeCSR->denseVec.begin() + control.rowsPerNode, nodeCSR->denseVec.end());
			} else if (control.myId < control.clusterCols && control.myId != 0) {
				// total number of rows in matrix not process or column
				MPI_Recv(&control.rowCount, 1, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
				usleep(100000 * control.myId);
				//std::cout << "rowCount = " << control.rowCount << std::endl;
				// Get rows and nnz per proc data
				nodeCSR->processData.resize(control.clusterRows*2);
				//std::cout << "receiving " << control.clusterRows*3 << " processData elements" << std::endl;
				MPI_Recv(&nodeCSR->processData[0], control.clusterRows*2, MPI_INT, 0, 0, control.row_comm,
				         MPI_STATUS_IGNORE);
				//std::cout << "processData.size() = " << nodeCSR->processData.size() << std::endl;

				//for (int i = 0; i < nodeCSR->processData.size(); i++){
				//	std::cout << nodeCSR->processData[i] << ", ";
				//}
				//std::cout << std::endl;

				control.elementCount = 0;
				for (int i = 0; i < control.clusterRows*2; i = i+2){
					control.elementCount += nodeCSR->processData[i];
				}
				//std::cout << "elementCount = " << control.elementCount << std::endl;
				int assignedRowCount = 0;
				for (int i = 1; i < control.clusterRows*2; i = i+2){
					assignedRowCount += nodeCSR->processData[i];
				}
				//std::cout << "assignedRowCount = " << assignedRowCount << std::endl;

				nodeCSR->assignedRowIds.resize(assignedRowCount);
				nodeCSR->csrRows.resize(assignedRowCount);
				nodeCSR->csrCols.resize(control.elementCount);
				nodeCSR->csrData.resize(control.elementCount);
				nodeCSR->denseVec.resize(control.rowCount);
				//std::cout << "Rows recieved: " << nodeCSR->csrRows.size() << ", NNZs received: " << nodeCSR->csrData.size() << std::endl;

				MPI_Recv(&nodeCSR->assignedRowIds[0], assignedRowCount, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrRows[0], assignedRowCount, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrCols[0], control.elementCount, MPI_INT, 0, 0, control.row_comm,
				         MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrData[0], control.elementCount, MPI_DOUBLE, 0, 0, control.row_comm,
				         MPI_STATUS_IGNORE);
				//std::cout << "Rows recieved: " << nodeCSR->csrRows.size() << ", NNZs received: " << nodeCSR->csrData.size() << ", denseVec received: " << nodeCSR->denseVec.size() << std::endl;
			}
			if (control.debug && control.myId == 0) std::cout << "Sending to column masters complete" << std::endl;
/*
			usleep(10000000 * control.myId);
			if (control.myId / control.clusterCols == 0){
				std::cout << "Rows recieved: " << nodeCSR->csrRows.size() << ", NNZs received: "
				          << nodeCSR->csrData.size() << ", denseVec received: " << nodeCSR->denseVec.size()
				          << std::endl;

				std::cout << "myId: " << control.myId << " - ";
				for (int i = 0; i < nodeCSR->processData.size(); i++){
					std::cout << nodeCSR->processData[i] << ", ";
				}
				std::cout << std::endl;
			}
*/

			// column masters send data to row nodes
			if (control.barrier) MPI_Barrier(control.col_comm);
			if (control.myId < control.clusterCols) {
				if (control.debug && control.myId == 0)	std::cout << "Column Masters Sending Data to col members" << std::endl;
				int rowsSent = 0, nnzSent = 0;
				for (int i = 1;
				     i < control.clusterRows; i++) {  // start at 1 since column master is the 0th node in the column
					rowsSent += nodeCSR->processData[((i-1)*2)+1];
					nnzSent += nodeCSR->processData[((i-1)*2)];

					MPI_Send(&control.rowCount, 1, MPI_INT, i, 0, control.col_comm);
					//std::cout << "sending " << nodeCSR->processData.size() << " processData elements to " << i << std::endl;
					MPI_Send(&(nodeCSR->processData[(i*2)]), 2, MPI_INT, i, 0, control.col_comm);
					//std::cout <<  "sent " << nodeCSR->processData[i*2] << ", " << nodeCSR->processData[(i*2)+1] << std::endl;

					MPI_Send(&(nodeCSR->assignedRowIds[rowsSent]), nodeCSR->processData[(i*2)+1], MPI_INT, i, 0,
					         control.col_comm);
					//std::cout << "sending csrRows to " << i << std::endl;
					MPI_Send(&(nodeCSR->csrRows[rowsSent]), nodeCSR->processData[(i*2)+1], MPI_INT, i, 0,
					         control.col_comm);
					//std::cout << "sending csrCols to " << i << std::endl;
					MPI_Send(&(nodeCSR->csrCols[nnzSent]), nodeCSR->processData[(i*2)], MPI_INT, i, 0,
					         control.col_comm);
					//std::cout << "sending csrData to " << i << std::endl;
					MPI_Send(&(nodeCSR->csrData[nnzSent]), nodeCSR->processData[(i*2)], MPI_DOUBLE, i, 0,
					         control.col_comm);
				}

				// Erase the excess data on the column master that has already been distributed to its row nodes
				int myLastData = nodeCSR->csrRows[control.rowsPerNode];
				if (!(nodeCSR->assignedRowIds.empty())) {
					nodeCSR->assignedRowIds.erase(nodeCSR->assignedRowIds.begin() + nodeCSR->processData[1], nodeCSR->assignedRowIds.end());
				}
				if (!(nodeCSR->csrRows.empty())) {
					nodeCSR->csrRows.erase(nodeCSR->csrRows.begin() + nodeCSR->processData[1], nodeCSR->csrRows.end());
				}
				if (!(nodeCSR->csrCols.empty())) {
					nodeCSR->csrCols.erase(nodeCSR->csrCols.begin() + nodeCSR->processData[0], nodeCSR->csrCols.end());
				}
				if (!(nodeCSR->csrData.empty())) {
					nodeCSR->csrData.erase(nodeCSR->csrData.begin() + nodeCSR->processData[0], nodeCSR->csrData.end());
				}
				//nodeCSR->rebase(control.myCol * control.colsPerNode);
			} else if (control.myId >= control.clusterCols) {
				nodeCSR->processData.resize(2,0);
				MPI_Recv(&control.rowCount, 1, MPI_INT, 0, 0, control.col_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&(nodeCSR->processData[0]), 2, MPI_INT, 0, 0, control.col_comm, MPI_STATUS_IGNORE);

				//usleep(10000000 * control.myId);
				//std::cout << "myId: " << control.myId << " - ";
				//for (int i = 0; i < nodeCSR->processData.size(); i++){
				//	std::cout << nodeCSR->processData[i] << ", ";
				//}
				//std::cout << std::endl;

				nodeCSR->assignedRowIds.resize(nodeCSR->processData[1]);
				nodeCSR->csrRows.resize(nodeCSR->processData[1]);
				nodeCSR->csrCols.resize(nodeCSR->processData[0]);
				nodeCSR->csrData.resize(nodeCSR->processData[0]);
				nodeCSR->denseVec.resize(control.rowCount);

				MPI_Recv(&nodeCSR->assignedRowIds[0], nodeCSR->processData[1], MPI_INT, 0, 0, control.col_comm,
				         MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrRows[0], nodeCSR->processData[1], MPI_INT, 0, 0, control.col_comm,
				         MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrCols[0], nodeCSR->processData[0], MPI_INT, 0, 0, control.col_comm,
				         MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrData[0], nodeCSR->processData[0], MPI_DOUBLE, 0, 0, control.col_comm,
				         MPI_STATUS_IGNORE);
				nodeCSR->rebase_balanced();
			}
			//std::cout << "Rows recieved: " << nodeCSR->csrRows.size() << ", NNZs received: " << nodeCSR->csrData.size() << ", denseVec received: " << nodeCSR->denseVec.size() << std::endl;
			if (control.myId == 0) nodeCSR->denseVec.resize(control.rowCount, 1.0);
			MPI_Bcast(&nodeCSR->denseVec[0], control.rowCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		dataTransmissionEnd = MPI_Wtime();

		//usleep(1000000 * control.myId);
		//std::cout << "myId: " << control.myId << " - ";
		//for (int i = 0; i < nodeCSR->processData.size(); i++){
		//	std::cout << nodeCSR->processData[i] << ", ";
		//}
		//std::cout << std::endl;

/*
			std::cout << "Process " << control.myId << std::endl;
			std::cout << "Rows: ";
			for (int j = 0; j < nodeCSR->csrRows.size(); j++){
				std::cout << nodeCSR->csrRows[j] << ",";
			}
			std::cout << std::endl << std::endl;

			std::cout << "Cols ";
			for (int j = 0; j < nodeCSR->csrCols.size(); j++){
				std::cout << nodeCSR->csrCols[j] << ",";
			}
			std::cout << std::endl << std::endl;

			std::cout << "Datas: ";
			for (int j = 0; j < nodeCSR->csrData.size(); j++){
				std::cout << nodeCSR->csrData[j] << ",";
			}
			std::cout << std::endl << std::endl;
*/

		// must be total number of rows since we want each process to take part in a collective reduction
		result.resize(control.rowCount, 0.0);

		if (control.debug && control.myId == 0) std::cout << "Starting SpMV computation" << std::endl;
		if (control.barrier) MPI_Barrier(MPI_COMM_WORLD);
		spmvStartTime = MPI_Wtime();
		if (nodeCSR->csrData.size() > 0) {
/*
			int rowErrorCount = 0, colErrorCount = 0;
			for (int i = 0; i < nodeCSR->csrRows.size(); i++) {
				if (nodeCSR->csrRows[i] > nodeCSR->csrCols.size()) {
					//std::cout << nodeCSR->csrRows[i] << std::endl;
					rowErrorCount++;
				}
			}
			for (int i = 0; i < nodeCSR->csrCols.size(); i++) {
				if (nodeCSR->csrCols[i] > control.rowCount) {
					//std::cout << "ERROR: INVALID DATA" << std::endl;
					colErrorCount++;
				}
			}
			if (rowErrorCount) std::cout << rowErrorCount << " row errors" << std::endl;
			if (colErrorCount) std::cout << colErrorCount << " col errors" << std::endl;
*/

			int ompThreadId, ompCPUId, start, end, i, j, k, rowsPerThread, rowEnd;
#pragma omp parallel num_threads(control.ompThreads) shared(nodeCSR, result) private(ompThreadId, ompCPUId, start, end, i, j, k, rowsPerThread, rowEnd)
			{
				ompThreadId = omp_get_thread_num();
				if (control.debug) {
					ompCPUId = sched_getcpu();
					//usleep(100000 * ompThreadId);
					std::cout << "Rank " << control.myId << ", Thread " << ompThreadId << " on core " << ompCPUId
					          << std::endl;
				}

				rowsPerThread = ceil(nodeCSR->csrRows.size() / (double) control.ompThreads);
				if (ompThreadId == control.ompThreads - 1) {
					rowEnd = nodeCSR->csrRows.size();
				} else {
					rowEnd = (ompThreadId + 1) * rowsPerThread;
				}

				if (ompThreadId == control.ompThreads - 1) {
					for (i = ompThreadId * rowsPerThread; i < nodeCSR->csrRows.size(); i++) {
						if (i == nodeCSR->csrRows.size() - 1) {
							for (j = nodeCSR->csrRows[i]; j < nodeCSR->csrData.size(); j++) {
								result[nodeCSR->assignedRowIds[i]] += nodeCSR->csrData[j] * nodeCSR->denseVec[nodeCSR->csrCols[j]];
							}
						} else {
							for (j = nodeCSR->csrRows[i]; j < nodeCSR->csrRows[i + 1]; j++) {
								result[nodeCSR->assignedRowIds[i]] += nodeCSR->csrData[j] * nodeCSR->denseVec[nodeCSR->csrCols[j]];
							}
						}
					}
				} else {
					for (i = ompThreadId * rowsPerThread; i < rowEnd; i++) {
						for (j = nodeCSR->csrRows[i]; j < nodeCSR->csrRows[i + 1]; j++) {
							result[nodeCSR->assignedRowIds[i]] += nodeCSR->csrData[j] * nodeCSR->denseVec[nodeCSR->csrCols[j]];
						}
					}
				}
			}
		}

		if (control.barrier) MPI_Barrier(MPI_COMM_WORLD);
		spmvEndTime = MPI_Wtime();
		if (control.debug && control.myId == 0) std::cout << "SpMV computation complete" << std::endl;

		if (control.masterOnly != true) {
			/*
			 *      MPI REDUCE FROM ALL TO MASTER
			 */
			if (control.debug && control.myId == 0) std::cout << "Starting MPI Reduction" << std::endl;
			reductionStartTime = MPI_Wtime();

			if (control.myId == 0) {
				MPI_Reduce(MPI_IN_PLACE, &result[0], control.rowCount, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			} else {
				MPI_Reduce(&result[0], &result[0], control.rowCount, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			}
			if (control.barrier) MPI_Barrier(MPI_COMM_WORLD);
			if (control.debug && control.myId == 0) std::cout << "MPI Reduction complete" << std::endl;
			reductionEndTime = MPI_Wtime();
		}

		if (control.myId == 0) {
			for (int i = 0; i < control.rowCount; i++){
				//std::cout << "result[" << i << "] = " << result[i] << std::endl;
			}
		}
	} else if (control.distributionMethod == 3){    // balanced distribution
		if (control.masterOnly != true) {
			if (control.debug && control.myId == 0) std::cout << "Sending data to column masters" << std::endl;
			// master to send data to cluster column masters
            //usleep(100000000 * control.myId);
			if (control.myId == 0) {
				for (int i = 1; i < control.clusterCols; i++) {  // start at 1 since Master is the row master
					MPI_Send(&control.rowCount, 1, MPI_INT, i, 0, control.row_comm);
					//MPI_Send(&control.maxRowsAssigned, 1, MPI_INT, i, 0, control.row_comm);
					//std::cout << "preparing to send " << clusterColData[i]->processData.size() << " processData elements to " << i << ": ";
					//for (int j = 0; j < clusterColData[i]->processData.size(); j++){
					//	std::cout << "for  = " << j << ": ";
					//	std::cout << clusterColData[i]->processData[j] << std::endl;
					//}
					//std::cout << std::endl;
					//std::cout  << "sending processData to " << i << std::endl;
					MPI_Send(&(clusterColData[i]->processData[0]), control.clusterRows*2, MPI_INT, i, 0, control.row_comm);
					//MPI_Send(&(clusterColData[i]->assignedRowIds[0]), clusterColData[i]->csrRows.size(), MPI_INT, i, 0, control.row_comm);
					//std::cout << "sending csrRows to " << i << std::endl;
					MPI_Send(&(clusterColData[i]->csrRows[0]), clusterColData[i]->csrRows.size(), MPI_INT, i, 0,
					         control.row_comm);
					//std::cout << "sending " << clusterColData[i]->csrCols.size() << " csrCols to " << i << std::endl;
					MPI_Send(&(clusterColData[i]->csrCols[0]), clusterColData[i]->csrCols.size(), MPI_INT, i, 0,
					         control.row_comm);
					//std::cout << "sending " << clusterColData[i]->csrData.size() << " csrData to " << i << std::endl;
					MPI_Send(&(clusterColData[i]->csrData[0]), clusterColData[i]->csrData.size(), MPI_DOUBLE, i, 0,
					         control.row_comm);
				}
				// delete and free column data that the master has already sent to the column masters, as it no longer needs
				// to be kept on the master
				for (int i = 1; i < clusterColData.size(); i++) {
					//delete (clusterColData[i]);
				}
				//clusterColData.erase(clusterColData.begin() + 1, clusterColData.end()); //remove transmitted data
				//nodeCSR->denseVec.erase(nodeCSR->denseVec.begin() + control.rowsPerNode, nodeCSR->denseVec.end());
			} else if (control.myId < control.clusterCols && control.myId != 0) {
				// total number of rows in matrix not process or column
				MPI_Recv(&control.rowCount, 1, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
				//MPI_Recv(&control.maxRowsAssigned, 1, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
				//std::cout << "rowCount = " << control.rowCount << std::endl;
				// Get rows and nnz per proc data
				nodeCSR->processData.resize(control.clusterRows*2);
				//std::cout << "receiving " << control.clusterRows*2 << " processData elements" << std::endl;
				MPI_Recv(&nodeCSR->processData[0], control.clusterRows*2, MPI_INT, 0, 0, control.row_comm,
				         MPI_STATUS_IGNORE);
				//std::cout << "processData.size() = " << nodeCSR->processData.size() << std::endl;

				//for (int i = 0; i < nodeCSR->processData.size(); i++){
				//	std::cout << nodeCSR->processData[i] << ", ";
				//}
				//std::cout << std::endl;

				control.elementCount = 0;
				for (int i = 0; i < control.clusterRows*2; i = i+2){
					control.elementCount += nodeCSR->processData[i];
				}
				//std::cout << "elementCount = " << control.elementCount << std::endl;
				int assignedRowCount = 0;
				for (int i = 1; i < control.clusterRows*2; i = i+2){
					assignedRowCount += nodeCSR->processData[i];
				}
				//std::cout << "assignedRowCount = " << assignedRowCount << std::endl;

				//nodeCSR->assignedRowIds.resize(assignedRowCount);
				nodeCSR->csrRows.resize(assignedRowCount);
				nodeCSR->csrCols.resize(control.elementCount);
				nodeCSR->csrData.resize(control.elementCount);
				nodeCSR->denseVec.resize(control.rowCount);
				//std::cout << "Rows recieved: " << nodeCSR->csrRows.size() << ", NNZs received: " << nodeCSR->csrData.size() << std::endl;

				//MPI_Recv(&nodeCSR->assignedRowIds[0], assignedRowCount, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrRows[0], assignedRowCount, MPI_INT, 0, 0, control.row_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrCols[0], control.elementCount, MPI_INT, 0, 0, control.row_comm,
				         MPI_STATUS_IGNORE);
				MPI_Recv(&nodeCSR->csrData[0], control.elementCount, MPI_DOUBLE, 0, 0, control.row_comm,
				         MPI_STATUS_IGNORE);
				//std::cout << "Rows recieved: " << nodeCSR->csrRows.size() << ", NNZs received: " << nodeCSR->csrData.size() << ", denseVec received: " << nodeCSR->denseVec.size() << std::endl;
			}

			if (control.myId == 0){
				nodeCSR->processData[0] = displacements[1];
				nodeCSR->processData[1] = rowCounts[0];
			}
			if (control.debug && control.myId == 0) std::cout << "Sending to column masters complete" << std::endl;
/*
			usleep(10000000 * control.myId);
			if (control.myId / control.clusterCols == 0){
				std::cout << "Rows recieved: " << nodeCSR->csrRows.size() << ", NNZs received: "
				          << nodeCSR->csrData.size() << ", denseVec received: " << nodeCSR->denseVec.size()
				          << std::endl;

				std::cout << "myId: " << control.myId << " - ";
				for (int i = 0; i < nodeCSR->processData.size(); i++){
					std::cout << nodeCSR->processData[i] << ", ";
				}
				std::cout << std::endl;
			}
*/

			// column masters send data to row nodes
			if (control.barrier) MPI_Barrier(control.col_comm);
			if (control.myId < control.clusterCols) {
				if (control.debug && control.myId == 0)	std::cout << "Column Masters Sending Data to col members" << std::endl;
				int rowsSent = 0, nnzSent = 0;
				for (int i = 1;
				     i < control.clusterRows; i++) {  // start at 1 since column master is the 0th node in the column
					rowsSent += nodeCSR->processData[((i-1)*2)+1];
					nnzSent += nodeCSR->processData[((i-1)*2)];

					MPI_Send(&control.rowCount, 1, MPI_INT, i, 0, control.col_comm);
					//std::cout << "sending maxRowsAssigned = " << control.maxRowsAssigned << " to " << i << std::endl;
					//MPI_Send(&control.maxRowsAssigned, 1, MPI_INT, i, 0, control.col_comm);
					//std::cout <<  "sending processdaata:  " << nodeCSR->processData[i*2] << ", " << nodeCSR->processData[(i*2)+1] << std::endl;
					MPI_Send(&(nodeCSR->processData[(i*2)]), 2, MPI_INT, i, 0, control.col_comm);
					//std::cout <<  "sent " << nodeCSR->processData[i*2] << ", " << nodeCSR->processData[(i*2)+1] << std::endl;
					//std::cout << "sending assigned Row Ids to " << i << std::endl;
					//MPI_Send(&(nodeCSR->assignedRowIds[rowsSent]), nodeCSR->processData[(i*2)+1], MPI_INT, i, 0,
					//         control.col_comm);
					//std::cout << "sending csrRows to " << i << std::endl;
					MPI_Send(&(nodeCSR->csrRows[rowsSent]), nodeCSR->processData[(i*2)+1], MPI_INT, i, 0,
					         control.col_comm);
					//std::cout << "sending csrCols to " << i << std::endl;
					MPI_Send(&(nodeCSR->csrCols[nnzSent]), nodeCSR->processData[(i*2)], MPI_INT, i, 0,
					         control.col_comm);
					//std::cout << "sending " << nodeCSR->processData[(i*2)] << " csrData to " << (control.myId*control.clusterCols)+(i*control.clusterRows) << std::endl;
					MPI_Send(&(nodeCSR->csrData[nnzSent]), nodeCSR->processData[(i*2)], MPI_DOUBLE, i, 0,
					         control.col_comm);
				}

				// Erase the excess data on the column master that has already been distributed to its row nodes
				int myLastData = nodeCSR->csrRows[control.rowsPerNode];
				if (!(nodeCSR->assignedRowIds.empty())) {
					nodeCSR->assignedRowIds.erase(nodeCSR->assignedRowIds.begin() + nodeCSR->processData[1], nodeCSR->assignedRowIds.end());
				}
				if (!(nodeCSR->csrRows.empty())) {
					nodeCSR->csrRows.erase(nodeCSR->csrRows.begin() + nodeCSR->processData[1], nodeCSR->csrRows.end());
				}
				if (!(nodeCSR->csrCols.empty())) {
					nodeCSR->csrCols.erase(nodeCSR->csrCols.begin() + nodeCSR->processData[0], nodeCSR->csrCols.end());
				}
				if (!(nodeCSR->csrData.empty())) {
					nodeCSR->csrData.erase(nodeCSR->csrData.begin() + nodeCSR->processData[0], nodeCSR->csrData.end());
				}
				//nodeCSR->rebase(control.myCol * control.colsPerNode);
			} else if (control.myId >= control.clusterCols) {
				nodeCSR->processData.resize(2,0);
				MPI_Recv(&control.rowCount, 1, MPI_INT, 0, 0, control.col_comm, MPI_STATUS_IGNORE);
				//std::cout << "received rowCount: " << control.rowCount << std::endl;
				//MPI_Recv(&control.maxRowsAssigned, 1, MPI_INT, 0, 0, control.col_comm, MPI_STATUS_IGNORE);
				//std::cout << "received maxRowsAssigned: " << control.maxRowsAssigned << std::endl;
				MPI_Recv(&(nodeCSR->processData[0]), 2, MPI_INT, 0, 0, control.col_comm, MPI_STATUS_IGNORE);
				//std::cout << "Proc " << control.myId << " received processData: " << nodeCSR->processData[0] << ", " << nodeCSR->processData[1] << std::endl;

				//nodeCSR->assignedRowIds.resize(nodeCSR->processData[1]);
				nodeCSR->csrRows.resize(nodeCSR->processData[1]);
				nodeCSR->csrCols.resize(nodeCSR->processData[0]);
				nodeCSR->csrData.resize(nodeCSR->processData[0]);
				nodeCSR->denseVec.resize(control.rowCount);

				//std::cout << "receiving assigned row ids" << std::endl;
				//MPI_Recv(&nodeCSR->assignedRowIds[0], nodeCSR->processData[1], MPI_INT, 0, 0, control.col_comm,
				//         MPI_STATUS_IGNORE);
				//std::cout << "receiving csrRows" << std::endl;
				MPI_Recv(&nodeCSR->csrRows[0], nodeCSR->processData[1], MPI_INT, 0, 0, control.col_comm,
				         MPI_STATUS_IGNORE);
				//std::cout << "receiving csrCols" << std::endl;
				MPI_Recv(&nodeCSR->csrCols[0], nodeCSR->processData[0], MPI_INT, 0, 0, control.col_comm,
				         MPI_STATUS_IGNORE);
				//std::cout << "receiving csrData" << std::endl;
				MPI_Recv(&nodeCSR->csrData[0], nodeCSR->processData[0], MPI_DOUBLE, 0, 0, control.col_comm,
				         MPI_STATUS_IGNORE);
                if(control.debug) std::cout << "Starting CSR rebasing" << std::endl;
				nodeCSR->rebase_balanced();
                if(control.debug) std::cout << "Done rebasing" << std::endl;
			}
			//std::cout << "Rows recieved: " << nodeCSR->csrRows.size() << ", NNZs received: " << nodeCSR->csrData.size() << ", denseVec received: " << nodeCSR->denseVec.size() << std::endl;
			if (control.myId == 0) nodeCSR->denseVec.resize(control.rowCount, 1.0);
            if(control.debug) std::cout << "Starting to broadcast dense vector" << std::endl;
			MPI_Bcast(&nodeCSR->denseVec[0], control.rowCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if(control.debug) std::cout << "Done broadcasting dense vector" << std::endl;
		}
		dataTransmissionEnd = MPI_Wtime();

		//usleep(1000000 * control.myId);
		//std::cout << "myId: " << control.myId << " - ";
		//for (int i = 0; i < nodeCSR->processData.size(); i++){
		//	std::cout << nodeCSR->processData[i] << ", ";
		//}
		//std::cout << std::endl;

/*
			std::cout << "Process " << control.myId << std::endl;
			std::cout << "Rows: ";
			for (int j = 0; j < nodeCSR->csrRows.size(); j++){
				std::cout << nodeCSR->csrRows[j] << ",";
			}
			std::cout << std::endl << std::endl;

			std::cout << "Cols ";
			for (int j = 0; j < nodeCSR->csrCols.size(); j++){
				std::cout << nodeCSR->csrCols[j] << ",";
			}
			std::cout << std::endl << std::endl;

			std::cout << "Datas: ";
			for (int j = 0; j < nodeCSR->csrData.size(); j++){
				std::cout << nodeCSR->csrData[j] << ",";
			}
			std::cout << std::endl << std::endl;
*/

		// must be total number of rows since we want each process to take part in a collective reduction
		if (control.myId == 0){
			result.resize(control.rowCount, 0.0);
		}
		if (control.myId != 0){
			gatheredResult.resize(nodeCSR->processData[1], 0.0);
		}
		std::cout << "gatheredResult.size() = " << gatheredResult.size() << std::endl;

		if (control.debug && control.myId == 0) std::cout << "Starting SpMV computation" << std::endl;
		if (control.barrier) MPI_Barrier(MPI_COMM_WORLD);
		spmvStartTime = MPI_Wtime();
		if (nodeCSR->csrData.size() > 0) {
/*
			int rowErrorCount = 0, colErrorCount = 0;
			for (int i = 0; i < nodeCSR->csrRows.size(); i++) {
				if (nodeCSR->csrRows[i] > nodeCSR->csrCols.size()) {
					//std::cout << nodeCSR->csrRows[i] << std::endl;
					rowErrorCount++;
				}
			}
			for (int i = 0; i < nodeCSR->csrCols.size(); i++) {
				if (nodeCSR->csrCols[i] > control.rowCount) {
					//std::cout << "ERROR: INVALID DATA" << std::endl;
					colErrorCount++;
				}
			}
			if (rowErrorCount) std::cout << rowErrorCount << " row errors" << std::endl;
			if (colErrorCount) std::cout << colErrorCount << " col errors" << std::endl;
*/

			int ompThreadId, ompCPUId, start, end, i, j, k, rowsPerThread, rowEnd;
			int totalRows = nodeCSR->csrRows.size();
			int totalData = nodeCSR->csrData.size();
			int dataStart, dataEnd;
#pragma omp parallel num_threads(control.ompThreads) shared(nodeCSR, gatheredResult, totalRows, totalData) private(ompThreadId, ompCPUId, start, end, i, j, k, rowsPerThread, rowEnd, dataStart, dataEnd)
			{
				ompThreadId = omp_get_thread_num();
				if (control.debug) {
					ompCPUId = sched_getcpu();
					//usleep(10000000 * control.myId );
					std::cout << "Rank " << control.myId << ", Thread " << ompThreadId << " on core " << ompCPUId
					          << std::endl;
				}

				rowsPerThread = ceil(nodeCSR->csrRows.size() / (double) control.ompThreads);
				if (ompThreadId == control.ompThreads - 1) {
					rowEnd = nodeCSR->csrRows.size();
				} else {
					rowEnd = (ompThreadId + 1) * rowsPerThread;
				}

				dataStart = nodeCSR->csrRows[i];
				dataEnd = nodeCSR->csrRows[i + 1];

				if (ompThreadId == control.ompThreads - 1) {
					for (i = ompThreadId * rowsPerThread; i < totalRows; i++) {
						if (i == totalRows - 1) {
							for (j = dataStart; j < totalData; j++) {
								gatheredResult[i] += nodeCSR->csrData[j] * nodeCSR->denseVec[nodeCSR->csrCols[j]];
							}
						} else {
							for (j = dataStart; j < dataEnd; j++) {
								gatheredResult[i] += nodeCSR->csrData[j] * nodeCSR->denseVec[nodeCSR->csrCols[j]];
							}
						}
					}
				} else {
					for (i = ompThreadId * rowsPerThread; i < rowEnd; i++) {
						for (j = dataStart; j < dataEnd; j++) {
							gatheredResult[i] += nodeCSR->csrData[j] * nodeCSR->denseVec[nodeCSR->csrCols[j]];
						}
					}
				}
			}
		}

		if (control.barrier) MPI_Barrier(MPI_COMM_WORLD);
		spmvEndTime = MPI_Wtime();
		if (control.debug && control.myId == 0) std::cout << "SpMV computation complete" << std::endl;

		if (control.masterOnly != true) {
			/*
			 *      MPI GATHER FROM ALL TO MASTER
			 */
			masterGatherStart = MPI_Wtime();
			if (control.debug && control.myId == 0){
				std::cout << "Starting MPI Gather" << std::endl;
			}
			//std::cout << gatheredResult.size() << ", " << rowCounts.size() << ", " << displacements.size() << std::endl;
			if (control.myId == 0) {
				//std::cout << "master in place " << rowCounts[0] << " rows" << std::endl;
				MPI_Gatherv(MPI_IN_PLACE, rowCounts[0], MPI_DOUBLE, &gatheredResult[0], &rowCounts[0], &displacements[0],
				           MPI_DOUBLE, 0, MPI_COMM_WORLD);
				//std::cout << " master done with gatherv" << std::endl;

			} else {
				//std::cout << control.myId << "sending " << nodeCSR->csrRows.size() << std::endl;
				MPI_Gatherv(&gatheredResult[0], nodeCSR->csrRows.size(), MPI_DOUBLE, &gatheredResult[0], &rowCounts[0],
				            &displacements[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				//std::cout << control.myId << " dont with gatherv" << std::endl;
			}
		}
		//std::cout << control.myId << " DONE WITH GATHER PORTION" << std::endl;
		if (control.barrier) MPI_Barrier(MPI_COMM_WORLD);

		if (control.myId == 0){
			//std::cout << "gatheredResult.size() = " << gatheredResult.size() << std::endl;
			for (int i = 1; i < control.rowDistribution.size(); i++) {
				//std::cout << "result[" << control.rowDistribution[i] << "]" << " = " << result[control.rowDistribution[i]] << " + " << gatheredResult[i] << std::endl;
				result[control.rowDistribution[i]] += gatheredResult[i];
			}
		}
		if (control.debug && control.myId == 0) std::cout << "MPI Gather complete" << std::endl;
		masterGatherEnd = MPI_Wtime();
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
	            //std::cout << masterData.result[i] << " = " << result[i] << std::endl;
                if (std::abs(masterData.result[i] - result[i]) > 0.001) {
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
		    if (control.distributionMethod == 1) {
			    std::cout << distributionEndTime - distributionStartTime << ","
			              << dataTransmissionEnd - dataTransmissionStart << "," << spmvEndTime - spmvStartTime << ","
			              << reductionEndTime - reductionStartTime << "," << masterGatherEnd - masterGatherStart << ","
			              << overallEndTime - overallStartTime << std::endl;
		    } else if (control.distributionMethod == 3){
			    std::cout << distributionEndTime - distributionStartTime << ","
			              << dataTransmissionEnd - dataTransmissionStart << "," << spmvEndTime - spmvStartTime << ","
			              << masterGatherEnd - masterGatherStart << "," << overallEndTime - overallStartTime << std::endl;
		    } else if (control.distributionMethod == 2){
			    std::cout << distributionEndTime - distributionStartTime << ","
			              << dataTransmissionEnd - dataTransmissionStart << "," << spmvEndTime - spmvStartTime << ","
			              << reductionEndTime - reductionStartTime << "," << overallEndTime - overallStartTime << std::endl;
		    }
	    }

    }
	if (control.debug && control.myId == 0) std::cout << "Output complete" << std::endl;


	return 0;
}
