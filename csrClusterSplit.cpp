//
// Created by brianpage on 5/15/17.
//
#include "csrClusterSplit.h"


//
//
//  Distribute sparse matrix amongst cluster nodes via OVERFLOW Distribution
//
//
void csrClusterSplit_Overflow(char *matrixFile, bool colMajor, std::string distributionMethod,
                              std::vector<int>& origin_row, std::vector<int>& origin_col,
                              std::vector<double>& origin_data, std::vector<std::vector<int> >& temp_row,
                              std::vector<std::vector<int> >& temp_col, std::vector<std::vector<double> >& temp_data,
                              std::vector<int>& colMasterTemp_row, std::vector<int>& colMasterTemp_col,
                              std::vector<double>& colMasterTemp_data, int& rowCount, int& colCount, int& nonZeros,
                              int& colsPerNode, int clusterRows, int clusterCols) {


    // out put converted format
    for (int i = 0; i < rowCount; i++){
        std::cout << "Row " << i << ": ";
        if (i != rowCount-1){
            for (int j = origin_row[i]; j < origin_row[i+1]; j++){
                std::cout << origin_data[j] << ", ";
            }
        } else {
            for (int j = origin_row[i]; j < origin_data.size(); j++){
                std::cout << origin_data[j] << ", ";
            }
        }
        std::cout << std::endl;
    }

    int maxRowLength = 0, rowLength = 0;
    for (int i = 0; i < rowCount; i++){
        int start = origin_row[i], stop = 0;

        if (i == rowCount-1){
            stop = origin_data.size();
        } else {
            stop = origin_row[i+1];
        }

        rowLength = stop-start;

        if (rowLength > maxRowLength){
            maxRowLength = rowLength;
        }
    }
    std::cout << "Max Row Length = " << maxRowLength << std::endl;

    colsPerNode = maxRowLength / clusterCols;
    if (colsPerNode < 1){
        colsPerNode = 1;
    }

    for (int i = 0; i < clusterCols; i++) {
        std::vector<int> a, b;
        std::vector<double> c;
        temp_row.push_back(a);
        temp_col.push_back(b);
        temp_data.push_back(c);
    }

    //std::cout << "Debug 1" << std::endl;
    for (int i = 0; i < rowCount; i++) {
        int start = origin_row[i];

        int stop;
        if (i == rowCount -1){
            stop = origin_data.size();
        } else {
            stop = origin_row[i+1];
        }

        int rowLength = stop - start;
        //std::cout  << "start = " << start << ", stop = " << stop << std::endl;
        int clusterColsNeededForRow = (rowLength/colsPerNode);
        for (int j = 0; j < clusterCols; j++) {
            if (j > clusterColsNeededForRow || rowLength <= j*colsPerNode) {
                //std::cout << "temp[" << j << "].push_back(-1)" << std::endl;
                temp_row[j].push_back(-1);
            } else {
                //std::cout << "temp[" << j << "].push_back(" << temp_data[j].size() << ")" << std::endl;
                temp_row[j].push_back(temp_data[j].size());
            }
        }
        int count = 0;
        for (int k = start; k < stop; k++){
            int col = count/colsPerNode;
            if (col == clusterCols){
                col--;
            }
            //std::cout << "adding " << origin_col[k] << " to temp_col[" << count/colsPerNode << "]" << std::endl;
            temp_col[col].push_back(origin_col[k]);
            //std::cout << "adding " << origin_data[k] << " to temp_data[" << count/colsPerNode << "]" << std::endl;
            temp_data[col].push_back(origin_data[k]);
            count++;
        }
    }
    std::cout << "DONE SPLITTING STUFF" << std::endl;
}


//
//
//  Distribute sparse matrix amongst cluster nodes via SPLIT MATRIX Distribution
//
//
void csrClusterSplit_SplitMatrix(controlData controlData, std::vector<csrSpMV*>& clusterColData){

    int colsLastColumn, colsPerColumn;
/*
    // out put converted format
    for (int i = 0; i < rowCount; i++){
        std::cout << "Row " << i << ": ";
        if (i != rowCount-1){
            for (int j = origin_row[i]; j < origin_row[i+1]; j++){
                std::cout << origin_data[j] << ", ";
            }
        } else {
            for (int j = origin_row[i]; j < origin_data.size(); j++){
                std::cout << origin_data[j] << ", ";
            }
        }
        std::cout << std::endl;
    }
*/

    clusterColData.resize(controlData.clusterCols);
    for (int i = 0; i < controlData.clusterCols; i++) {
        clusterColData[i] = new csrSpMV;
    }

    //
    // Read MMF file and insert elements into clusterColData as needed
    //
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

                if (!(tempRow == 0 && tempCol == 0 && tempData == 0.0)) {
                    elements.push_back(tempElement);
                }
            }
        } else {
            previousLineCommented = true;
        }
    }
    printf("%d rows, %d cols, and %d non-zeros\n", rowCount, colCount, nonZeros);

    // sort the vector of elements by row, and then each row, based on column
    //std::stable_sort(elements.begin(), elements.end(), sortByCol);

    // add to csr vectors for use as CSR Format
    for (int k = 0; k < elements.size(); k++) {
        //std::cout << "k = " << k << std::endl;

        if (k == 0){
            //std::cout << k << ", " << elements[k].col << ", " << elements[k].data << std::endl;
            csr_row.push_back(k);
        } else {
            if (elements[k].col != elements[k-1].col){
                //std::cout << k << ", " << elements[k].col << ", " << csr_row[rowsAdded-1] << std::endl;

                csr_row.push_back(k);
            }
        }

        csr_col[k] = elements[k].row;
        csr_data[k] = elements[k].data;
    }

    /*
    for (int k = 0; k < csr_row.size(); k++) {
        std::cout << csr_row[k] << std::endl;
    }

    for (int k = 0; k < csr_data.size(); k++) {
        std::cout << csr_data[k] << std::endl;
    }
    */

    if (rowCount != csr_row.size()){
        std::cout << "Actual row count does NOT match reported rowCount" << std::endl;
    }

    //
    // DONE
    //


















    colsPerColumn = controlData.rowCount / controlData.clusterCols;
    if (controlData.rowCount % controlData.clusterCols != 0){
        colsLastColumn = colsPerColumn + (controlData.rowCount % controlData.clusterCols);
    } else {
        colsLastColumn = controlData.rowCount / controlData.clusterCols;
    }


    for (int i = 0; i < rowCount; i++) {
        int start, stop;

        start = origin_row[i];
        if (i == rowCount-1){
            stop = origin_data.size();
        } else {
            stop = origin_row[i+1];
        }

        //std::cout << "start = " << start << ", stop = " << stop << std::endl;

        int rowStarts[clusterCols];
        for (int c = 0; c < clusterCols; c++){
            rowStarts[c] = -1;
        }

        for (int j = start; j < stop; j++){
            //indicate that the row actually has something in it for this column of cluster nodes
            //std::cout << "origin_col[" << j << "] = " << origin_col[j] << std::endl;
            int col = origin_col[j]/(rowCount / clusterCols);
            if (col == clusterCols){
                col -= 1;
            }
            if (rowStarts[col] == -1){
                //std::cout << "rowStarts[" << col << "] = " << colMasterTemp_data[col].size() << std::endl;
                rowStarts[col] = colMasterTemp_data[col].size();
            }

            //std::cout << "colMasterTemp_col[" << col << "].push_back(" << origin_col[j] << ")" << std::endl;
            colMasterTemp_col[col].push_back(origin_col[j]);
            //std::cout << "colMasterTemp_data[" << col << "].push_back(" << origin_data[j] << ")" << std::endl;
            colMasterTemp_data[col].push_back(origin_data[j]);
        }

        for (int j = 0; j < clusterCols; j++){
            //std::cout << "colMasterTemp_row[" << j << "].push_back(" << rowStarts[j] << ")" << std::endl;
            colMasterTemp_row[j].push_back(rowStarts[j]);
            //std::cout << "colMasterTemp_row[" << j << "][" << i << "] = " << colMasterTemp_row[j][colMasterTemp_row[j].size()-1] << std::endl;
            //std::cout << "colMasterTemp_row[" << j << "].size() = " << colMasterTemp_row[j].size() << std::endl;
        }

    }

    std::cout << "DONE SPLITTING STUFF" << std::endl;
}




//
//
//  Distribute sparse matrix amongst cluster nodes via BALANCED Distribution
//
//
void csrClusterSplit_ElementBalanced(char *matrixFile, bool colMajor, std::string distributionMethod, int processCount,
                                     std::vector<int>& origin_row, std::vector<int>& origin_col,
                                     std::vector<double>& origin_data, std::vector<std::vector<int> >& temp_row,
                                     std::vector<std::vector<int> >& temp_col,
                                     std::vector<std::vector<double> >& temp_data, int& rowCount, int& colCount,
                                     int& nonZeros, int& colsPerNode, int clusterRows,
                                     int clusterCols, std::vector<std::vector<std::vector <int> > >& nodeRowOwnership){

    std::cout << "Inside csrClusterSplit_ElementBalanced()" << std::endl;


    /*
    // out put converted format
    for (int i = 0; i < rowCount; i++){
        std::cout << "Row " << i << ": ";
        if (i != rowCount-1){
            for (int j = origin_row[i]; j < origin_row[i+1]; j++){
                std::cout << origin_data[j] << ", ";
            }
        } else {
            for (int j = origin_row[i]; j < origin_data.size(); j++){
                std::cout << origin_data[j] << ", ";
            }
        }
        std::cout << std::endl;
    }
    */

    colsPerNode = rowCount / clusterRows;
    int rowOverFlow = rowCount % colsPerNode;
    int nodeBalanceElementCount = nonZeros / processCount;
    std::cout << "nonZeros = " << nonZeros << ", processCount = " << processCount << std::endl;
    std::cout << "Optimal elements per node = " << nodeBalanceElementCount << std::endl;

    std::vector<std::vector <int> > rowLengths;
    int elementCount = 0;
    //std::cout << "about to add rows to rowLengths vector" << std::endl;
    for (int i = 0; i < rowCount; i++){
        std::vector<int> temp(3);
        if(i == rowCount-1){
            temp[0] = origin_data.size() - origin_row[i];
        } else {
            temp[0] = origin_row[i + 1] - origin_row[i];
        }
        temp[1] = i;
        temp[2] = elementCount;

        //std::cout << "about to push_back temp onto rowLenghts" << std::endl;
        rowLengths.push_back(temp);
        elementCount += temp[0];
    }

    // find the balanced element distribution for the given cluster and matrix size
    balanceDistribution(processCount, nodeBalanceElementCount, rowLengths, nodeRowOwnership);

    // add rows and elements to the proper temp vectors
    for (int i = 0; i < processCount; i++) {
        std::vector<int> rows, cols;
        std::vector<double> data;

        for (int j = 0; j < nodeRowOwnership[i][1].size(); j++){    // for every row assigned to process
            if (data.size() == 0){
                rows.push_back(0);
            } else {
                rows.push_back(data.size());
            }

            int start, stop;
            //if (j == nodeRowOwnership[i][1].size()-1){
                start = nodeRowOwnership[i][3][j];
                stop = nodeRowOwnership[i][3][j] + nodeRowOwnership[i][2][j];
            //} else {
            //    start = nodeRowOwnership[i][3][j];
            //    stop = nodeRowOwnership[i][3][j+1];
            //}

            for (int k = start; k < stop; k++){
                //std::cout << "i = " << i << ", j = " << j << ", k = " << k << ", start = " << start << ", stop = " << stop
                //          << std::endl;
                cols.push_back(origin_col[k]);
                data.push_back(origin_data[k]);
            }
        }

        temp_row.push_back(rows);
        temp_col.push_back(cols);
        temp_data.push_back(data);
    }
}



void balanceDistribution(int processCount, int nodeBalanceElementCount, std::vector<std::vector <int> >& rowLengths,
                         std::vector<std::vector <std::vector <int> > >& tempPacking){

    //std::cout << "Inside balanceDitribution()" << std::endl;
    //std::cout << "rowlenghts.size() = " << rowLengths.size() << std::endl;

    //
    // do a bin packing/set partitioning of the rows up to nodeBalanceElementCount. For the bin packing, put the rows
    // into nodeRowOwnership and update the first element of that vector with the sum of the elements given to that node
    //

    // If any row exceeds the balanced element count, split the row into smaller chunks. That is, we create new "rows"
    // such that their contents are the separate chunks of the original row.
    for (int i = 0; i < rowLengths.size(); i++){
        if (rowLengths[i][0] > nodeBalanceElementCount){
            std::vector <int> temp(3);
            int splitCount =  rowLengths.size()/nodeBalanceElementCount;
            int rowStart = rowLengths[i][2];

            if (splitCount > 1) { // more than one split needed
                for (int j = 0; j < splitCount; j++) {
                    if (j == splitCount-1){
                        temp[0] = rowLengths[i][0] % nodeBalanceElementCount; // last split gets overflow
                    } else {
                        temp[0] = nodeBalanceElementCount;
                    }
                    temp[1] = rowLengths[i][1];
                    temp[2] = rowStart +(j*nodeBalanceElementCount);
                    rowLengths.push_back(temp);
                }
            } else {
                temp[0] = rowLengths[i][0] - nodeBalanceElementCount;   // give split row the overflow
                temp[1] = rowLengths[i][1]; // make sure new split row has the same row number
                temp[2] = rowStart + nodeBalanceElementCount;
                rowLengths.push_back(temp);
            }
        }
    }


    std::sort(rowLengths.begin(), rowLengths.end());    // sort by the length of each row or split row

    for (int i = 0; i < rowLengths.size()/2; i++){
        std::vector <int> temp = rowLengths[i];

        rowLengths[i] = rowLengths[rowLengths.size()-(i+1)];
        rowLengths[rowLengths.size()-(i+1)] = temp;
    }

    /*
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "----- RowLengths has the following rows/data -----" << std::endl;
    for (int i = 0; i < rowLengths.size(); i++){
        std::cout << "Row " << rowLengths[i][1] << " - Length = " << rowLengths[i][0] << " starting at " << rowLengths[i][2] << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    */

    // we now have rows no larger than nodeBalanceElement's value and the rows have been sorted according to their
    // lengths. Next we want to perform a bin packing in an effort to balance the number of elements as evenly as
    // possible across all nodes in the cluster

    /*
     *  tempPacking and nodeRowOwnership are meant to have the following structure:
     *
     *      tempPacking[i][j][k]
     *
     *                  i -> 0 to p and represents cluster nodes based on MPI process ID
     *
     *                      j -> 0 to 3: 0 = sum of elements assigned to ith node
     *                                   1 = vector containing the row IDs owned by the ith node
     *                                   2 = vector containing element count for each row in element of j=1
     *                                   3 = vector containing starting element id of row or row split in to origin_data
     *
     *                          k -> 0 to however many rows assigned to the ith node based on balancing algorithm
     *
     */
    //initialize tempPacking
    for (int i = 0; i < processCount; i++){
        std::vector<int> temp1dVec;
        std::vector<std::vector <int> > temp2dVec;

        for(int j = 0; j < 4; j++) {
            if (j == 0) {
                std::vector<int> startingAssignedElements;
                startingAssignedElements.push_back(0);
                temp2dVec.push_back(startingAssignedElements);
            } else {
                temp2dVec.push_back(temp1dVec);
            }
        }
        tempPacking.push_back(temp2dVec);
    }

    //std::cout << 1 << std::endl;
/*

    // Greedy set partitioning using AVL tree to keep track of total elements assigned to each node as we assign rows or
    // row splits to a node, we update that nodes tree node with the new element count
    nodeTree elementToClusterTree;

    // add nodes for all MPI processes (cluster nodes) to tree
    for (int i = 1; i < processCount; i++) {
        std::cout << "adding " << i << " to tree" << std::endl;
        elementToClusterTree.insert(0, i);
    }

    //elementToClusterTree.balance(elementToClusterTree.root);
    std::cout << "Tree Height = " << elementToClusterTree.height(elementToClusterTree.root) << std::endl;
    elementToClusterTree.inorder(elementToClusterTree.root);
    elementToClusterTree.display(elementToClusterTree.root,0);
    std::cout << std::endl;

    elementToClusterTree.display(elementToClusterTree.root,1);
    std::cout << std::endl;

    int nodeAssignment;
    //std::cout << 3 << std::endl;
    for (int i = 0; i < rowLengths.size(); i++){
        //std::cout << 4 << std::endl;
        nodeAssignment = -1;
        //std::cout << 5 << std::endl;
        elementToClusterTree.assignElements(elementToClusterTree.root, rowLengths[i][0], nodeAssignment);
        //std::cout << 6 << std::endl;

        if (nodeAssignment == -1){
            std::cout << "ERROR: Node assignment error during element distribution balancing" << std::endl;
            exit(1);
        }
        //std::cout << 7 << std::endl;

        // update tempPacking with new row assignment to cluster node
        tempPacking[nodeAssignment][0][0] += rowLengths[i][0];
        tempPacking[nodeAssignment][1].push_back(rowLengths[i][1]);
        tempPacking[nodeAssignment][2].push_back(rowLengths[i][0]);
        tempPacking[nodeAssignment][3].push_back(rowLengths[i][2]);
        //std::cout << 8 << std::endl;
    }
    //std::cout << 9 << std::endl;
    elementToClusterTree.freeTree();
*/

    int nodeAssignment;
    //std::cout << 3 << std::endl;
    for (int i = 0; i < rowLengths.size(); i++){
        //std::cout << "i = " << i << std::endl;

        int nodeAssignment = 0, min = tempPacking[nodeAssignment][0][0];
        for (int j = 0; j < processCount; j++){
            if (tempPacking[j][0][0] < min){
                min = tempPacking[j][0][0];
                nodeAssignment = j;
            }
        }
        //std::cout << 6 << std::endl;

        // update tempPacking with new row assignment to cluster node
        tempPacking[nodeAssignment][0][0] += rowLengths[i][0];
        tempPacking[nodeAssignment][1].push_back(rowLengths[i][1]);
        tempPacking[nodeAssignment][2].push_back(rowLengths[i][0]);
        tempPacking[nodeAssignment][3].push_back(rowLengths[i][2]);
        //std::cout << 8 << std::endl;
    }

    /*
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "----- TempPacking Layout -----" << std::endl;
    for (int i = 0; i < tempPacking.size(); i++){
        std::cout << "Process " << i << std::endl;
        for (int j = 0; j < tempPacking[i][1].size(); j++){
            std::cout << "Row " << tempPacking[i][1][j] << ", length " << tempPacking[i][2][j] << " starts at " << tempPacking[i][3][j] << std::endl;
        }
        std::cout << std::endl;

    }
    std::cout << std::endl;
    std::cout << std::endl;
    */
}

void splitDenseVector_ElementBalanced(std::vector<double> denseVector,
                                      std::vector<std::vector<double> >& splitDenseVector,
                                      std::vector<std::vector<std::vector <int> > >& nodeRowOwnership){

    for (int i = 0; i < nodeRowOwnership.size(); i++){
        std::vector<double> temp;
        for (int j = 0; j < nodeRowOwnership[i][1].size(); j++){
            temp.push_back(denseVector[j]);
        }
        splitDenseVector.push_back(temp);
    }

}
