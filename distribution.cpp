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
		//std::cout << control.rowCount << "," << control.colCount << "," << control.nonZeros << std::endl;
		if (line[0] != '%') {
			//std::cout << line << std::endl;
			if (previousLineCommented == true) {
				//some stuff to get row,col and size data, etc.
				previousLineCommented = false;

				size_t pos = 0;
				std::string token;
				i = 0;
				while ((pos = line.find(' ')) != std::string::npos) {
					token = line.substr(0, pos);
					line.erase(0, pos + 1);

					//std::cout << "token = " << token << std::endl;

					if (i == 0) {
						control.rowCount = std::stoi(token);
					} else {
						control.colCount = std::stoi(token);
					}

					i++;
				}
				//std::cout << "token = " << line << std::endl;
				control.nonZeros = std::stoi(line);

				control.rowsPerNode = ceil(control.rowCount / (float) control.clusterRows);
				control.colsPerNode = ceil(control.colCount / (float) control.clusterCols);

				std::cout << "rowsPerNode = " << control.rowsPerNode << "colsPerNode = " << control.colsPerNode << std::endl;

				//std::cout << control.rowCount << "," << control.colCount << "," << control.nonZeros << std::endl;

				// determine where to test overflow for extra rows that must be added to the last cluster row's work
				if (control.rowCount % control.clusterCols != 0) {
					rowsLastRow = control.rowsPerNode + (control.rowCount % control.clusterRows);
				} else {
					rowsLastRow = control.rowCount / control.clusterCols;
				}
				control.lastClusterRowRowStart = control.colCount - rowsLastRow;
				// determine where to test overflow for extra columns that must be added to the last cluster column's
				// work
				if (control.colCount % control.clusterCols != 0) {
					colsLastColumn = control.colsPerNode + (control.colCount % control.clusterCols);
				} else {
					colsLastColumn = control.colCount / control.clusterCols;
				}
				control.lastClusterColColStart = control.colCount - colsLastColumn;

				//std::cout << "lastClusterColColStart = " << control.lastClusterColColStart << std::endl;

			} else {
				//Element tempElement(0, 0, 0.0);
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
				if (!(tempRow < 0 || tempCol < 0 || tempData == 0.0)) {
					elements.emplace_back(tempRow, tempCol, tempData);
				}

				//std::cout << "assignedCol = " << assignedCol << std::endl;
				// if the row is
				//
				/*
				if (tempRow == previousRow) {
					// Do nothing, row is already present in list
				} else {
					// Add new row to cluster Column, as it has not been seen previously
					for (int k = 0; k < clusterColData.size(); k++) {
						clusterColData[k]->csrRows.push_back(clusterColData[k]->csrData.size());
					}
				}


				// assign the column and data to the correct csrSpMV object for the cluster column it belongs to
				clusterColData[assignedCol]->csrCols.push_back(tempCol);
				clusterColData[assignedCol]->csrData.push_back(tempData);

				previousRow = tempRow;

				//std::cout << tempRow << "," << tempCol << "," << tempData << std::endl;
				*/
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

		// if the number of columns does not evenly divide amongst the number of cluster columns, the final
		// cluster column is given the excess whereas all other columns receive the same amount of columns to
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


    //std::cout << "out of pushing onto distribtuon" << std::endl;

    for (int i = 0; i < control.clusterCols; i ++){
        if (clusterColData[i]->csrRows.size() != control.rowCount){
            clusterColData[i]->csrRows.resize(control.rowCount, (clusterColData[i]->csrData.size()-1));
        }
        //std::cout << "clusterColData[i]->csrRows.size() = " << clusterColData[i]->csrRows.size() << std::endl;
    }

    //std::cout << "out of resizing csrRows" << std::endl;

    for (int i = 0; i < control.rowCount; i++){
        //clusterColData[0]->denseVec.push_back((double) (rand()) / (double) (RAND_MAX));
	    clusterColData[0]->denseVec.push_back(1.0);
    }

    //std::cout << "out of pushing onto denseVec" << std::endl;

    // Output the number of Non Zero elements that have been assigned to each cluster column, as well as each cluster
	// node within that column. This will be used for memory transfer and performance evaluation calculations
/*
	std::cout << std::endl << "Distribution of NonZero Elements" << std::endl;
    	for (int i = 0; i < control.clusterCols; i++){
        	std::cout << "Column " << i << ": " << clusterColData[i]->csrData.size() << std::endl;

	    	for (int j = 0; j < control.clusterRows; j++) {
                std::cout << "\tRow " << j << ": ";

                int start = clusterColData[i]->csrRows[j * control.rowsPerNode];
                int end;

                if (j == control.clusterRows - 1) {
                    end = clusterColData[i]->csrData.size();
                } else {
                    end = clusterColData[i]->csrRows[(j + 1) * control.rowsPerNode];
                }
                std::cout << end - start << std::endl;
                std::cout << "start = " << start << ", end = " << end << std::endl;
	    	}
    	}
*/

/*
	std::cout << std::endl << "Distribution of NonZero Elements" << std::endl;
    	for (int i = 0; i < control.clusterCols; i++){
        	std::cout << "Column " << i << ": " << clusterColData[i]->csrData.size() << std::endl;
            for (int j = 0; j < clusterColData[i]->csrData.size(); j++){
                std::cout << clusterColData[i]->csrData[j] << ",";
            }
            std::cout << std::endl;
   	}
*/
	
 }
