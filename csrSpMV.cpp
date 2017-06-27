//
// Created by brianpage on 6/22/17.
//

#include "csrSpMV.h"

struct Element {
    int row;
    int col;
    double data;

    Element(int r, int c, double d) : row(r), col(c), data(d) {}
};

csrSpMV::~csrSpMV() {
    csrRows.clear();
    csrCols.clear();
    csrData.clear();
    denseVec.clear();
    result.clear();
}

void csrSpMV::nodeSpMV(controlData controlData1) {


}

void csrSpMV::masterOnlySpMV(controlData controlData) {
    //convert sparse matrix from Matrix Market to Compressed Sparse Row format
    std::vector <Element> elements; //holds each matrix element as read from the Matrix Market format file
    int tempRow, tempCol;
    double tempData;

    // Read in sparse matrix saved in Matrix Market Format
    std::ifstream infile(controlData.matrixFile);
    if (!infile){
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
                        controlData.rowCount = std::stoi(token);
                    } else {
                        controlData.colCount = std::stoi(token);
                    }

                    i++;
                }

                controlData.nonZeros = std::stoi(line);
            } else {
                Element tempElement(0, 0, 0.0);
                size_t pos = 0;
                std::string token;
                i = 0;
                while ((pos = line.find(' ')) != std::string::npos) {
                    token = line.substr(0, pos);
                    line.erase(0, pos + 1);

                    if (i == 0) {
                        tempElement.row = ::atoi(token.c_str()) - 1;
                    } else {
                        tempElement.col = ::atoi(token.c_str()) - 1;
                    }

                    i++;
                }
                tempElement.data = ::atof(line.c_str());

                if (!(tempElement.row == 0 && tempCol == 0 && tempData == 0.0)) {
                    elements.push_back(tempElement);
                }
            }
        } else {
            previousLineCommented = true;
        }
    }
    printf("%d rows, %d cols, and %d non-zeros\n", controlData.rowCount, controlData.colCount, controlData.nonZeros);

    // add to csr vectors for use as CSR Format
    csrCols.resize(controlData.nonZeros);
    csrData.resize(controlData.nonZeros);

    int rowsAdded = -1;
    for (int k = 0; k < elements.size(); k++) {
        //std::cout << "k = " << k << std::endl;

        if (k == 0){
            //std::cout << k << ", " << elements[k].col << ", " << elements[k].data << std::endl;
            csrRows.push_back(k);
            rowsAdded++;
        } else {
            if (elements[k].col != elements[k-1].col){
                //std::cout << k << ", " << elements[k].col << ", " << csr_row[rowsAdded-1] << std::endl;

                csrRows.push_back(k);
                rowsAdded++;
            }
        }

        csrCols[k] = elements[k].row;
        csrData[k] = elements[k].data;
        //std::cout << "Added " << csr_data[k] << " to row " << csr_row.size()-1 << std::endl;
    }

    if (controlData.rowCount != csrRows.size()){
        std::cout << "Actual row count does NOT match reported rowCount" << std::endl;
    }


    //std::cout << "Creating Dense Vector" << std::endl;
    if (!controlData.myId) {
        //std::cout << "Populating Dense Vector" << std::endl;
        for (int i = 0; i < controlData.rowCount; i++) {
            denseVec.push_back(1.0);
        }
    }

    std::cout << "Performing Master Only SpMV" << std::endl;
    result.resize(controlData.rowCount);

    std::cout << "rowcount = " << controlData.rowCount << ", csrRows.size() = " << csrRows.size() << std::endl;
    for (int i = 0; i < controlData.rowCount; i++) {
        //std::cout << "i = " << i << std::endl;
        int temp = 0.0;
        if (i != controlData.rowCount - 1) {
            if (controlData.colMajor) { //col major order selected
                for (int j = csrRows[i]; j < csrRows[i + 1]; j++) {   // go to end of current row
                    // entire row is multiplied by a single dense vector element
                    //std::cout << "i = " << i << " j = " << j << ", origin_row["<< i << "] = " << origin_row[i] << std::endl;
                    temp += csrData[j] * denseVec[i];
                }
            } else {    // row major order selected
                for (int j = csrRows[i]; j < csrRows[i + 1]; j++) {   // go to end of current row
                    //std::cout << "i = " << i << " j = " << j << ", origin_row["<< i << "] = " << origin_row[i] << std::endl;
                    temp += csrData[j] * denseVec[csrCols[j]];
                }
            }
        } else {
            if (controlData.colMajor) { //col major order selected
                for (int j = csrRows[i]; j < csrData.size(); j++) {  // go to end of data vector
                    // entire row is multiplied by a single dense vector element
                    //std::cout << "i = " << i << " j = " << j << ", origin_row["<< i << "] = " << origin_row[i] << std::endl;
                    temp += csrData[j] * denseVec[i];
                }
            } else {    // row major order selected
                for (int j = csrRows[i]; j < csrData.size(); j++) {  // go to end of data vector
                    //std::cout << "i = " << i << " j = " << j << ", origin_row["<< i << "] = " << origin_row[i] << std::endl;
                    temp += csrData[j] * denseVec[csrCols[j]];
                }
            }
        }
        result[i] = temp;
    }
}

