//
// Created by brianpage on 5/15/17.
//
#include "csrClusterSplit.h"

void csrClusterSplit_Overflow(char *matrixFile, bool colMajor, std::string distributionMethod,
                              std::vector<int>& origin_row, std::vector<int>& origin_col,
                              std::vector<double>& origin_data, std::vector<std::vector<int> >& temp_row,
                              std::vector<std::vector<int> >& temp_col, std::vector<std::vector<double> >& temp_data,
                              std::vector<int>& colMasterTemp_row, std::vector<int>& colMasterTemp_col,
                              std::vector<double>& colMasterTemp_data, int& rowCount, int& colCount, int& nonZeros,
                              int& colsPerNode, int clusterRows, int clusterCols) {

    //convert sparse matrix from Matrix Market to Compressed Sparse Row format
    if (colMajor) {
        std::cout << "*** Using column major order of input matrix ***" << std::endl;
        MMCOO_to_CSR_colMajor(matrixFile, origin_row, origin_col, origin_data, rowCount, colCount, nonZeros);
    } else {
        MMCOO_to_CSR(matrixFile, origin_row, origin_col, origin_data, rowCount, colCount, nonZeros);
    }

    colsPerNode = rowCount / clusterRows;
    int rowOverFlow = rowCount % colsPerNode;

    for (int i = 0; i < clusterCols; i++) {
        std::vector<int> a, b;
        std::vector<double> c;
        temp_row.push_back(a);
        temp_col.push_back(b);
        temp_data.push_back(c);
    }

    for (int i = 0; i < rowCount; i++) {
        int firstData = origin_row[i];
        int lastData, rowElements;

        if (i == rowCount - 1) {
            lastData = origin_data.size() - 1;
            rowElements = origin_data.size() - firstData;
        } else {
            lastData = origin_row[i + 1] - 1;
            rowElements = lastData - firstData;
        }

        int maxMatCols = colsPerNode;

        int clusterColsNeeded = rowElements / maxMatCols;
        for (int k = 0; k < clusterCols; k++) {
            if (k <= clusterColsNeeded) {
                temp_row[k].push_back(temp_data[k].size());
            } else {
                temp_row[k].push_back(-1);
            }
        }

        int count = 0;
        for (int j = firstData; j <= lastData; j++) {
            temp_data[count / maxMatCols].push_back(origin_data[j]);
            temp_col[count / maxMatCols].push_back(origin_col[j]);
            count++;
        }
    }

    for (int k = 0; k < temp_data[0].size(); k++) {
        colMasterTemp_data.push_back(temp_data[0][k]);
        colMasterTemp_col.push_back(temp_col[0][k]);
    }
    for (int k = 0; k < rowCount; k++) {
        colMasterTemp_row.push_back(temp_row[0][k]);
    }


    //std::cout << "rowCount = " << rowCount << ", nonZeros = " << nonZeros  << std::endl;
    //std::cout << "origin_row.size() = " << origin_row.size() << ", origin_col.size() = " << origin_col.size() << ", origin_data.size() = " << origin_data.size() << std::endl;
    //std::cout << "temp_row.size() = " << temp_row.size() << ", temp_col.size() = " << temp_col.size() << ", temp_data.size() = " << temp_data.size() << std::endl;
}

void csrClusterSplit_ElementBalanced(char *filename, bool colMajor, std::string distributionMethod, int processCount,
                                     std::vector<int>& origin_row, std::vector<int>& origin_col,
                                     std::vector<double>& origin_data, std::vector<std::vector<int> >& temp_row,
                                     std::vector<std::vector<int> >& temp_col,
                                     std::vector<std::vector<double> >& temp_data, std::vector<int>& colMasterTemp_row,
                                     std::vector<int>& colMasterTemp_col, std::vector<double>& colMasterTemp_data,
                                     int& rowCount, int& colCount, int& nonZeros, int& colsPerNode, int clusterRows,
                                     int clusterCols, std::vector<std::vector <int> >& nodeRowOwnership){

    colsPerNode = rowCount / clusterRows;
    int rowOverFlow = rowCount % colsPerNode;
    int nodeBalanceElementCount = nonZeros / processCount;

    std::vector<std::vector <int> > rowLengths;
    for (int i = 0; i < rowCount; i++){
        std::vector<int> temp(2);
        temp[0] = origin_row.size();
        temp[1] = i;
    }

    std::sort(rowLengths.begin(), rowLengths.end());
    balanceDistribution(processCount, nodeBalanceElementCount, rowLengths, nodeRowOwnership);

    for (int i = 0; i < processCount; i++) {
        std::vector<int> a, b;
        std::vector<double> c;
        temp_row.push_back(a);
        temp_col.push_back(b);
        temp_data.push_back(c);
    }
    //
    // add rows and elements to the proper temp vectors here
    //
}

void balanceDistribution(int processCount, int nodeBalanceElementCount, std::vector<std::vector <int> >& nodeElementCount,
                         std::vector<std::vector <int> >& nodeRowOwnership){
    //
    // do a bin packing of the rows up to nodeBalanceElementCount. For the bin packing, put the rows into
    // nodeRowOwndership and update the first element of that vector with the sum of the elements given to that node
    //
}