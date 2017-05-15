//
// Created by brianpage on 5/15/17.
//
#include "clusterSpMV.h"

void clusterSpMV(int ompThreads, std::vector<int> csr_row, std::vector<int> csr_col, std::vector<double> csr_data,
                 double *denseVector, double *nodeResult, int rowsPerNode) {

    int i, j, rowStart, rowEnd, nextValidRow, ompId, firstRow, lastRow;
/*
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
*/

    #pragma omp parallel num_threads(ompThreads) shared(nodeResult, denseVector, csr_row, csr_data, csr_col) \
        private(ompId, firstRow, lastRow, rowStart, rowEnd, nextValidRow, i, j)
    {
        ompId = omp_get_thread_num();

        if (ompThreads < csr_row.size()) { // we want this
            std::cout << "outside ***************, ompThreads = " << ompThreads << std::endl;
            if (csr_row.size() % ompThreads == 0) {
                firstRow = ompId * csr_row.size() / ompThreads;
                lastRow = (ompId + 1) * csr_row.size() / ompThreads;
            } else {
                firstRow = ompId * csr_row.size() / ompThreads;
                if (ompId != ompThreads - 1) {
                    lastRow = (ompId + 1) * csr_row.size() / ompThreads;
                } else {
                    lastRow = csr_row.size();   // doesnt exist but we use less than this value anyways
                }
            }
        } else {
            std::cout << "inside, ompThreads = " << ompThreads  << std::endl;
             if (ompId < csr_row.size()) {
                firstRow = ompId;
                if (ompId == csr_row.size() - 1) {
                    lastRow = csr_row.size(); // doesnt exist but we use less than this value anyways
                } else {
                    lastRow = ompId + 1;
                }
            }
        }
        //std::cout << "Thread " << ompId << ", firstRow = " << firstRow << ", lastRow = " << lastRow << std::endl;

        if (ompId < rowsPerNode) {
            std::cout << "Thread " << ompId << ", firstRow = " << firstRow << ", lastRow = " << lastRow << std::endl;

            for (i = firstRow; i < lastRow; i++) { // iterate through all rows node is to work on
                if (csr_row[i] != -1) { // check is i points to row that has no data on this node
                    if (i == (rowsPerNode - 1)) { // check if i is last row of node
                        for (j = csr_row[i] - csr_row[0];
                             j < csr_data.size(); j++) { // go from last row start to end of data
                            //#pragma omp atomic
                            nodeResult[i] += csr_data[j] * denseVector[csr_col[j]];
                        }
                    } else {
                        nextValidRow = i + 1;
                        while ((csr_row[nextValidRow] == -1) && (nextValidRow < csr_row.size())) {
                            nextValidRow++;
                        }

                        rowStart = csr_row[i] - csr_row[0];
                        rowEnd = csr_row[nextValidRow] - csr_row[0];
                        for (int j = rowStart; j < rowEnd; j++) { // go from last row start to before row i+1
                            //#pragma omp atomic
                            nodeResult[i] += csr_data[j] * denseVector[csr_col[j]];
                        }
                    }
                }
            }
            //#pragma omp barrier
        } else {
            //std::cout << myId << "-" << ompId << " doing nothing" << std::endl;
        }
        //std::cout << ompId << " Done!" << std::endl;
    }   /* END of omp parallel for  */
}