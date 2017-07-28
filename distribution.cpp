//
// Created by brianpage on 5/15/17.
//
#include "distribution.h"

//
//
//  Distribute sparse matrix amongst cluster nodes via SPLIT MATRIX Distribution
//
//
void distribution_SplitMatrix(controlData& control, std::vector<csrSpMV*>& clusterColData) {
    int rowsLastRow, rowsPerRow, colsLastColumn, colsPerColumn;

    clusterColData.resize(control.clusterCols);
    for (int i = 0; i < control.clusterCols; i++) {
        clusterColData[i] = new csrSpMV;
    }

    //
    // Read MMF file and insert elements into clusterColData as needed
    //
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

                control.rowsPerNode = ceil(control.rowCount / (float)control.clusterRows);
                control.colsPerNode = ceil(control.colCount / (float)control.clusterCols);

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
            } else {
                size_t pos = 0;
                std::string token;
                i = 0;
                while ((pos = line.find(' ')) != std::string::npos) {
                    token = line.substr(0, pos);
                    line.erase(0, pos + 1);

                    if (i == 0) {
                        tempCol = ::atoi(token.c_str()) - 1;
                    } else {
                        tempRow = ::atoi(token.c_str()) - 1;
                    }

                    i++;
                }
                tempData = ::atof(line.c_str());

                // if the number of columns does not evenly divide amongst the number of cluster columns, the final
                // cluster column is given the excess whereas all other columns receive the same amount of columns to
                // work over.
                if (tempCol > control.lastClusterColColStart) {
                    assignedCol = control.clusterCols - 1;
                } else {
                    assignedCol = tempCol / control.colsPerNode;
                }

                //std::cout << tempCol << ", " << tempRow << ", " << tempData << ", assignedCol = " << assignedCol
                //fa          << std::endl;

                // if the row is
                //
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

            }
        } else {
            previousLineCommented = true;
        }
    }

    for (int i = 0; i < control.clusterCols; i ++){
        if (clusterColData[i]->csrRows.size() != control.rowCount){
            clusterColData[i]->csrRows.resize(control.rowCount, (clusterColData[i]->csrData.size()-1));
        }
    }

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
	    }
    }

	/*
	for (int i = 0; i < control.clusterCols; i++){
        std::cout << "ID = " << i << std::endl;
        for (int j = 0; j < clusterColData[i]->csrRows.size(); j++) {
            std::cout << "Row " << j << ": ";

            int start = clusterColData[i]->csrRows[j];
            int end;

            if (j == clusterColData[i]->csrRows.size() - 1) {
                end = clusterColData[i]->csrData.size();
            } else {
                end = clusterColData[i]->csrRows[j + 1];
            }

            for (int k = start; k < end; k++) {
                std::cout << clusterColData[i]->csrData[k] << ", ";
            }
            std::cout << std::endl;
        }

        std::cout << std::endl << std::endl;
    }
    */
 }
