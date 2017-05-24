//
// Created by brianpage on 4/4/17.
//
#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>

struct Element {
    int row;
    int col;
    double data;

    Element(int r, int c, double d) : row(r), col(c), data(d) {}
};

bool sortByCol(const Element& lhs, const Element& rhs) {
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

void MMCOO_to_CSR(char *filename, std::vector<int> &csr_row, std::vector<int> &csr_col,
                  std::vector<double> &csr_data, int &rowCount, int &colCount, int &nonZeros) {

    std::vector <Element> elements; //holds each matrix element as read from the Matrix Martket format file
    int tempRow, tempCol;
    double tempData;

    // Read in sparse matrix saved in Matrix Market Format
    std::ifstream infile(filename);
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
                        rowCount = std::stoi(token);
                    } else {
                        colCount = std::stoi(token);
                    }

                    i++;
                }

                nonZeros = std::stoi(line);
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
    printf("%d rows, %d cols, and %d non-zeros\n", rowCount, colCount, nonZeros);

    // sort the vector of elements by row, and then each row, based on column
    std::stable_sort(elements.begin(), elements.end(), sortByCol);

    // add to csr vectors for use as CSR Format
    int previousRow = -1;
    for (i = 0; i < nonZeros; i++) {
        csr_data.push_back(elements[i].data);
        csr_col.push_back(elements[i].col);

        if (previousRow == -1) {
            csr_row.push_back(elements[i].row);
        } else {
            if (previousRow == elements[i].row) {

            } else {
                csr_row.push_back(i);
            }
        }

        previousRow = elements[i].row;
    }
}

void MMCOO_to_CSR_colMajor(char *matrixFile, std::vector<int> &csr_row, std::vector<int> &csr_col,
                  std::vector<double> &csr_data, int &rowCount, int &colCount, int &nonZeros) {

    std::vector <Element> elements; //holds each matrix element as read from the Matrix Market format file
    int tempRow, tempCol;
    double tempData;

    // Read in sparse matrix saved in Matrix Market Format
    std::ifstream infile(matrixFile);
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
                        rowCount = std::stoi(token);
                    } else {
                        colCount = std::stoi(token);
                    }

                    i++;
                }

                nonZeros = std::stoi(line);
            } else {
                Element tempElement(0, 0, 0.0);
                size_t pos = 0;
                std::string token;
                i = 0;
                while ((pos = line.find(' ')) != std::string::npos) {
                    token = line.substr(0, pos);
                    line.erase(0, pos + 1);

                    if (i == 0) {
                        tempElement.col = ::atoi(token.c_str()) - 1;    // rows become columns
                    } else {
                        tempElement.row = ::atoi(token.c_str()) - 1;    // columns become rows
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
    printf("%d rows, %d cols, and %d non-zeros\n", rowCount, colCount, nonZeros);

    // sort the vector of elements by row, and then each row, based on column
    std::stable_sort(elements.begin(), elements.end(), sortByCol);

    // add to csr vectors for use as CSR Format
    int previousRow = -1;
    for (i = 0; i < nonZeros; i++) {
        csr_data.push_back(elements[i].data);
        csr_col.push_back(elements[i].col);

        if (previousRow == -1) {
            csr_row.push_back(elements[i].row);
        } else {
            if (previousRow == elements[i].row) {

            } else {
                csr_row.push_back(i);
            }
        }

        previousRow = elements[i].row;
    }

    if (rowCount != csr_row.size()){
        std::cout << "Actual row count does NOT match reported rowCount" << std::endl;
    }
}