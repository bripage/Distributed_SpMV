//
// Created by brianpage on 5/15/17.
//
#include "clusterSpMV.h"

void clusterSpMV_overflow(int ompThreads, std::vector<int> csr_row, std::vector<int> csr_col, std::vector<double> csr_data,
                 double *denseVector, double *nodeResult, int rowsPerNode, bool colMajor) {

    int i, j, rowStart, rowEnd, nextValidRow, ompId, firstRow, lastRow;
    for (i = 0; i < rowsPerNode; i++) { // iterate through all rows node is to work on
        if (csr_row[i] != -1) { // check is i points to row that has no data on this node
            if (i == (rowsPerNode) - 1) { // check if i is last row of node
                for (j = csr_row[i] - csr_row[0];
                     j < csr_col.size(); j++) { // go from last row start to end of data
                    nodeResult[i] += csr_data[j] * denseVector[csr_col[j]];
                }
            } else {
                nextValidRow = i + 1;
                while ((csr_row[nextValidRow] == -1) && (nextValidRow < csr_row.size())) {
                    nextValidRow++;
                }
                rowStart = csr_row[i] - csr_row[0];
                rowEnd = csr_row[nextValidRow] - csr_row[0];
                for (j = rowStart; j < rowEnd; j++) { // go from last row start to before row i+1
                    nodeResult[i] += csr_data[j] * denseVector[csr_col[j]];
                }
            }
        }
    }


}

void clusterSpMV_SplitMatrix(int ompThreads, std::vector<int> csr_row, std::vector<int> csr_col,
                             std::vector<double> csr_data, std::vector<double> denseVector,
                             std::vector<double>& nodeResult, bool colMajor){


}

void clusterSpMV_ElementBalanced(int ompThreads, std::vector<int> csr_row, std::vector<int> csr_col,
                                 std::vector<double> csr_data, std::vector<double> denseVector,
                                 std::vector<double>& nodeResult, bool colMajor){

    int i, j, rowStart, rowEnd, nextValidRow, ompId, firstRow, lastRow;

    for (i = 0; i < csr_row.size(); i++) { // iterate through all rows node is to work on
        rowStart = csr_row[i];
        if (i == csr_row.size()-1){
            rowEnd = csr_data.size();
        } else {
            rowEnd = csr_row[i+1];
        }

        //std::cout << "i = " << i << ", rowStart = " << rowStart << ", rowEnd = " << rowEnd << std::endl;
        for (j = rowStart; j < rowEnd; j++) {
            nodeResult[i] += csr_data[j] * denseVector[i];
            //std::cout << "nodeResult[" << i << "] += " << csr_data[j] << " * " << denseVector[i] << " --> nodeResult["
            //          << i << "] = " << nodeResult[i] << std::endl;
        }
    }

}