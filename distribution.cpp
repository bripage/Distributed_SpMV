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
		clusterColData[i]->processRowCounts.resize(control.clusterRows, 0);
	}


	std::vector <std::vector <Element> > distributionRows; //holds each matrix element as read from the Matrix Market format file
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
				for (int i = 0; i < control.rowCount; i++){
					distributionRows[i] = new std::vector<elements>;
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
					distributionRows[tempCol].emplace_back(tempRow, tempCol, tempData);
				}
			}
		} else {
			previousLineCommented = true;
		}
	}

	//
	// Greedy packing method to determine "balanced" NNZ distribution
	//
	int avgNNZperProcess = control.nonZeros / control.processCount;
	int maxRowSize = 0; // this will determine how much any single row can throw off the dirstributions balance
	for (int i = 0; i < control.rowCount; i++){
		if (maxRowSize < distributionRows[i].size()){
			maxRowSize = distributionRows[i].size();
		}
	}

	std::vector <std::vector <int>> procRowAssignment;
	for (int i = 0; i < control.processCount; i++){
		std::vector<int> temp;
		procRowAssignment.push_back(temp);
	}

	for (int i = 0; i < control.clusterCols; i++){
		for (int j = 0; j < control.clusterRows; j++) {
			std::vector<int> temp;

			bool filled = false;
			do {
				for (int k = 0; j < distributionRows.size(); k++) {
					if (distributionRows[k][0] != -1) {
						if (clusterColData[i]->processRowCounts[j] < avgNNZperProcess) {
							if ((avgNNZperProcess - clusterColData[i]->processRowCounts[j]) > ((clusterColData[i]->processRowCounts[j] + distributionRows[k].size())-avgNNZperProcess)) {
								procRowAssignment[(j * control.clusterCols) + i].push_back(k);
								clusterColData[i]->processRowCounts[j] += distributionRows[i].size();
							}
						} else {
							filled = true;
							break;
						}
					}
				}
			} while (!filled);

		}
	}


	for (int i = 0; i < control.processCount; i++){
		std::cout << clusterColData[i%control.clusterCols]->processRowCounts[i/control.clusterRows] << std::endl;
	}
}