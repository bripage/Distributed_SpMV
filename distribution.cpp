//
// Created by brianpage on 5/15/17.
//
#include "distribution.h"

bool sortRowThenCol(const Element& lhs, const Element& rhs) {
	if (lhs.row == rhs.row) {
		if (lhs.col > rhs.col) {
			return false;
		} else {
			return true;
		}
	} else if ( lhs.row > rhs.row){
		return false;
	} else {
		return true;
	}
}

bool sortByLength(const row &a, const row &b){
	return a.rowLength < b.rowLength;
}

//
//  Distribute sparse matrix amongst cluster nodes via SPLIT MATRIX Distribution
//
void distribution_SplitMatrix(controlData& control, std::vector<csrSpMV*>& clusterColData) {
	srand(static_cast <unsigned> (time(0)));
	int rowsLastRow, rowsPerRow, colsLastColumn, colsPerColumn;

	//
	// Read MMF file and insert elements into clusterColData as needed
	//
	clusterColData.resize(control.clusterCols);
	std::vector <Element> elements; //holds each matrix element as read from the Matrix Market format file
	std::vector <int> rowNNZs, colNNZs; // holds the number of NNZ per row for statistic calculation if enabled

	for (int i = 0; i < control.clusterCols; i++) {
		clusterColData[i] = new csrSpMV;
	}

	int tempRow, tempCol, assignedCol = 0;
	double tempData;

	// Read in sparse matrix saved in Matrix Market Format
	std::ifstream infile(control.matrixFile);
	if (!infile) {
		std::cout << "FAILED TO OPEN FILE!" << std::endl;
		exit(1);
	}
	std::string line;

	int i = 0, j = 0, previousRow = -1;
	bool previousLineCommented = false;
	while (std::getline(infile, line)) {
		if (line[0] != '%') {
			if (previousLineCommented == true) {
				//some stuff to get row,col and size data, etc.
				previousLineCommented = false;

				size_t pos = 0;
				std::string token;
				i = 0;
				while ((pos = line.find(' ')) != std::string::npos) {
					token = line.substr(0, pos);
					line.erase(0, pos + 1);

					if (i == 0) {
						control.rowCount = std::stoi(token);
					} else {
						control.colCount = std::stoi(token);
					}
					i++;
				}

				control.nonZeros = std::stoi(line);
				control.rowsPerNode = ceil(control.rowCount / (float) control.clusterRows);
				control.colsPerNode = ceil(control.colCount / (float) control.clusterCols);
				if (control.debug && control.myId == 0) std::cout << "RowCount = " << control.rowCount
				                                                  << ", ColCount = " << control.colCount
				                                                  << ", NNZ Count = " << control.nonZeros << std::endl;
				if (control.debug && control.myId == 0) std::cout << "RowsPerProcess = " << control.rowsPerNode
				                                                  << ", ColsPerProcess = " << control.colsPerNode
				                                                  << std::endl;

				// determine where to test overflow for extra rows that must be added to the last cluster row's work
				if (control.rowCount % control.clusterCols != 0) {
					rowsLastRow = control.rowsPerNode + (control.rowCount % control.clusterRows);
				} else {
					rowsLastRow = control.rowCount / control.clusterCols;
				}
				control.lastClusterRowRowStart = control.colCount - rowsLastRow;
				// determine where to test overflow for extra columns that must be added to the last cluster column's
				// work
				control.lastClusterColColStart = control.colsPerNode*(control.clusterCols-1);
			} else {
				size_t pos = 0;
				std::string token;
				i = 0;
				while ((pos = line.find(' ')) != std::string::npos) {
					token = line.substr(0, pos);
					line.erase(0, pos + 1);

					if (i == 0) {
						tempRow = ::atoi(token.c_str()) - 1;
					} else {
						tempCol = ::atoi(token.c_str()) - 1;
					}

					i++;
				}
				tempData = ::atof(line.c_str());

				// If we have read a valid element data create an element object for it
				if (!(tempData == 0.0 || tempData == 0)) {
					elements.emplace_back(tempRow, tempCol, tempData);
				}
			}
		} else {
			previousLineCommented = true;
		}
	}

	//
	//  Sort the elements by row and then by column
	//
	std::stable_sort(elements.begin(), elements.end(), sortRowThenCol);
	for (int i = 0; i < elements.size(); i++){
		//make sure all columns have a row pointer for the current elements row. If not create it!
		if (elements[i].row == previousRow) {
			// Do nothing, row is already present in list
		} else {
			// Add new row to cluster Column, as it has not been seen previously
			for (int k = 0; k < clusterColData.size(); k++) {
				clusterColData[k]->csrRows.push_back(clusterColData[k]->csrData.size());
			}
		}

		// if the number of columns does not evenly divide amongst the number of process columns, the final
		// process column is given the excess whereas all other columns receive the same amount of columns to
		// work over.
		if (elements[i].col > control.lastClusterColColStart) {
			assignedCol = control.clusterCols - 1;
		} else {
			assignedCol = elements[i].col / control.colsPerNode;
		}

		// assign the element to the correct process matrix column it belongs to
		clusterColData[assignedCol]->csrCols.push_back(elements[i].col);
		clusterColData[assignedCol]->csrData.push_back(elements[i].data);

		previousRow = elements[i].row;
	}

    for (int i = 0; i < control.clusterCols; i ++){
        if (clusterColData[i]->csrRows.size() != control.rowCount){
            clusterColData[i]->csrRows.resize(control.rowCount, (clusterColData[i]->csrData.size()-1));
        }
    }

    for (int i = 0; i < control.rowCount; i++){
        //clusterColData[0]->denseVec.push_back((double) (rand()) / (double) (RAND_MAX));
	    clusterColData[0]->denseVec.push_back(0.1234567);
    }

 }

//
//  Distribute sparse matrix amongst cluster nodes via Balanced NNZ per process Distribution
//
void distribution_Balanced(controlData& control, std::vector<csrSpMV*>& clusterColData) {
	//
	// Read MMF file and insert elements into clusterColData as needed
	//
	clusterColData.resize(control.clusterCols);
	for (int i = 0; i < control.clusterCols; i++) {
		clusterColData[i] = new csrSpMV;
		clusterColData[i]->processData.resize(control.clusterRows*3, 0);
	}


	std::vector <row> distributionRows; //holds each matrix element as read from the Matrix Market format file
	std::vector <int> rowNNZs, colNNZs; // holds the number of NNZ per row for statistic calculation if enabled

	int tempRow, tempCol, assignedCol = 0;
	double tempData;

	// Read in sparse matrix saved in Matrix Market Format
	std::ifstream infile(control.matrixFile);
	if (!infile) {
		std::cout << "FAILED TO OPEN FILE!" << std::endl;
		exit(1);
	}
	std::string line;

	int i = 0, j = 0, previousRow = -1;
	bool previousLineCommented = false;
	while (std::getline(infile, line)) {
		if (line[0] != '%') {
			if (previousLineCommented == true) {
				//some stuff to get row,col and size data, etc.
				previousLineCommented = false;

				size_t pos = 0;
				std::string token;
				i = 0;
				while ((pos = line.find(' ')) != std::string::npos) {
					token = line.substr(0, pos);
					line.erase(0, pos + 1);

					if (i == 0) {
						control.rowCount = std::stoi(token);
					} else {
						control.colCount = std::stoi(token);
					}
					i++;
				}

				control.nonZeros = std::stoi(line);
				if (control.debug && control.myId == 0) std::cout << "RowCount = " << control.rowCount
				                                                  << ", ColCount = " << control.colCount
				                                                  << ", NNZ Count = " << control.nonZeros << std::endl;
				if (control.debug && control.myId == 0) std::cout << "RowsPerProcess = " << control.rowsPerNode
				                                                  << ", ColsPerProcess = " << control.colsPerNode
				                                                  << std::endl;
				distributionRows.resize(control.rowCount);
				for (int i = 0; i < control.rowCount; i++){
					distributionRows[i].rowLength = 0;
					distributionRows[i].processAssignment = -1;
					distributionRows[i].rowId = i;
				}
			} else {
				size_t pos = 0;
				std::string token;
				i = 0;
				while ((pos = line.find(' ')) != std::string::npos) {
					token = line.substr(0, pos);
					line.erase(0, pos + 1);

					if (i == 0) {
						tempRow = ::atoi(token.c_str()) - 1;
					} else {
						tempCol = ::atoi(token.c_str()) - 1;
					}

					i++;
				}
				tempData = ::atof(line.c_str());

				// If we have read a valid element data create an element object for it
				if (!(tempData == 0.0 || tempData == 0)) {
					distributionRows[tempCol].rowIds.push_back(tempRow);
					distributionRows[tempCol].data.push_back(tempData);
					distributionRows[tempCol].rowLength++;
				}
			}
		} else {
			previousLineCommented = true;
		}
	}

	if (control.debug) std::cout << "Splitting rows of excess length" << std::endl;
	// split rows that are longer the avgNNZperProcess
	int avgNNZperProcess = ceil((double)control.nonZeros / (double)control.processCount);
	for (int i = 0; i < distributionRows.size(); i++){
		if (distributionRows[i].rowLength > avgNNZperProcess){
			row splitRow;

			//create and fill new row from large row split
			splitRow.rowId = distributionRows[i].rowId;
			for (int j = avgNNZperProcess; j < distributionRows[i].data.size(); j++){
				splitRow.rowIds.push_back(distributionRows[i].rowIds[j]);
				splitRow.data.push_back(distributionRows[i].data[j]);
			}
			splitRow.rowLength = splitRow.data.size();
			splitRow.processAssignment = -1;

			//remove copied elements from original row that was split by this process
			distributionRows[i].rowIds.erase(distributionRows[i].rowIds.begin() + avgNNZperProcess, distributionRows[i].rowIds.end());
			distributionRows[i].data.erase(distributionRows[i].data.begin() + avgNNZperProcess, distributionRows[i].data.end());
			distributionRows[i].rowLength = distributionRows[i].data.size();
		}
	}
	if (control.debug) std::cout << "Done spliting rows" << std::endl;

	if (control.debug) std::cout << "Sorting Rows" << std::endl;
	//sort rows based on row length
	if (control.debug) std::sort(distributionRows.begin(), distributionRows.end(), sortByLength);
	if (control.debug) std::cout << "Done sorting rows" << std::endl;

	if (control.debug) std::cout << "Assigning rows to processes" << std::endl;
	//
	// Greedy packing method to determine "balanced" NNZ distribution
	//
	std::vector <int> nnzAssignedPerProc(control.processCount, 0);
	for (int i = 0; i < control.processCount; i++){
		//std::cout << "i = " << i << std::endl;
		bool filled = false;
		while (!filled) {
			for (int j = 0; j < distributionRows.size(); j++) {
				if (distributionRows[j].processAssignment == -1) {
					if (nnzAssignedPerProc[i] < avgNNZperProcess) {
						//std::cout << "i = " << i << ", j = " << j << ", nnzAssignedPerProc[" << i << "] = " << nnzAssignedPerProc[i] << ", avgNNZperProcess = " << avgNNZperProcess << std::endl;
						if ((avgNNZperProcess - nnzAssignedPerProc[i]) > ((nnzAssignedPerProc[i] + distributionRows[j].rowLength) - avgNNZperProcess)) {
							distributionRows[j].processAssignment = i;
							nnzAssignedPerProc[i] += distributionRows[j].rowLength;
						}
					} else {
						filled = true;
						break;
					}
				}
				if (j == distributionRows.size() - 1){
					filled = true;
				}
			}
		}
	}
	if (control.debug) std::cout << "Done assigning rows" << std::endl;

	for (int i = 0; i < control.processCount; i++){
		std::cout << nnzAssignedPerProc[i] << std::endl;
	}

	if (control.debug) std::cout << "Populating clusterColData" << std::endl;
	for (int i = 0; i < control.processCount; i++){
		for (int j = 0; j < distributionRows.size(); j++){
			if (distributionRows[i].processAssignment == i){
				std::cout << " istributionRows[" << i << "].processAssignment = " << distributionRows[i].processAssignment << std::endl;
				clusterColData[i%control.clusterCols]->processData[((i/control.clusterRows)*3)+1]+=1;
				std::cout << "clusterColData[" << i%control.clusterCols << "]->processData["
				          << ((i/control.clusterRows)*3)+1 << "] = "
				          << clusterColData[i%control.clusterCols]->processData[((i/control.clusterRows)*3)+1]
				          << std::endl;
				clusterColData[i%control.clusterCols]->csrRows.push_back(clusterColData[i%control.clusterCols]->csrData.size());
				for (int k = 0; k < distributionRows[j].data.size(); k++) {
					clusterColData[i % control.clusterCols]->csrCols.push_back(distributionRows[j].rowIds[k]);
					clusterColData[i % control.clusterCols]->csrData.push_back(distributionRows[j].data[k]);
				}
				clusterColData[i%control.clusterCols]->processData[((i/control.clusterRows)*3)] += distributionRows[j].rowLength;
				//clusterColData[i%control.clusterCols]->denseVec.push_back(denseVec[j]);
				clusterColData[i%control.clusterCols]->denseVec.push_back(0.1234567);
				clusterColData[i%control.clusterCols]->processData[((i/control.clusterRows)*3)+2]+=1;
			}
		}
	}
	if (control.debug) std::cout << "Done populating clusterColData" << std::endl;
}