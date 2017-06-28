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

    int i = 0, j = 0;
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

                printf("%d rows, %d cols, and %d non-zeros\n", control.rowCount, control.colCount, control.nonZeros);

                control.rowsPerNode = ceil(control.rowCount / (float)control.clusterRows);
                control.colsPerNode = ceil(control.colCount / (float)control.clusterCols);
                std::cout << "rowsPerNode = " << control.rowsPerNode << ", colsPerNode = " << control.colsPerNode << std::endl;

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

                std::cout << " rowsLastRow = " << colsLastColumn << ", lastClusterRowRowStart = " << control.lastClusterRowRowStart << std::endl;
                std::cout << " colsLastColumn = " << colsLastColumn << ", lastClusterColColStart = " << control.lastClusterColColStart << std::endl;
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

                // if the number of columns does not event divide amongst the number of cluster columns, the final
                // cluster column is given the excess whereas all other columns receive the same amount of columns to
                // work over.
                if (tempCol > control.lastClusterColColStart) {
                    assignedCol = control.clusterCols - 1;
                } else {
                    assignedCol = tempCol / control.colsPerNode;
                }

                // if the row is
                //
                if (clusterColData[assignedCol]->csrRows.empty()) {
                    clusterColData[assignedCol]->csrRows.push_back(clusterColData[assignedCol]->csrData.size());
                } else {
                    if (tempRow == clusterColData[assignedCol]->csrRows.back()) {
                        // Do nothing, row is already present in list
                    } else {
                        // Add new row to cluster Column, as it has not been seen previously
                        clusterColData[assignedCol]->csrRows.push_back(clusterColData[assignedCol]->csrData.size());
                    }
                }

                // assign the column and data to the correct csrSpMV object for the cluster column it belongs to
                clusterColData[assignedCol]->csrCols.push_back(tempCol);
                clusterColData[assignedCol]->csrData.push_back(tempData);

            }
        } else {
            previousLineCommented = true;
        }
    }

    std::cout << "DONE SPLITTING MATRIX AMONGST CLUSTER COLUMNS" << std::endl;


    for (int i = 0; i < control.clusterCols; i++){
        std::cout << "Column " << i << ": " << clusterColData[i]->csrData.size() << std::endl;
    }
}
