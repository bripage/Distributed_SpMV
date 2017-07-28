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
csrSpMV::csrSpMV(){

}

csrSpMV::csrSpMV(const csrSpMV& objToCopy){
	csrRows = objToCopy.csrRows;
	csrCols = objToCopy.csrCols;
	csrData = objToCopy.csrData;
}


csrSpMV::~csrSpMV() {
    csrRows.clear();
    csrCols.clear();
    csrData.clear();
    denseVec.clear();
    result.clear();
}

void csrSpMV::nodeSpMV(controlData control, std::vector <double>& result) {

    //std::this_thread::sleep_for(std::chrono::milliseconds(control.myId * 2000));
    //std::cout << "myID = " << control.myId << ", data length = " << csrData.size() << std::endl;

    if (csrData.size() > 0) {
		int ompThreadId, start, end, i, j;

	    std::cout << "omp_get_max_threads() = " << omp_get_max_threads() << std::endl;

		#pragma omp parallel num_threads(control.ompThreads) shared(result) private(ompThreadId, start, end, i, j)
	    {
		    ompThreadId = omp_get_thread_num();

		    std::cout << "Thread " << ompThreadId << " starting " << std::endl;

	        for (i = 0; i < csrRows.size(); i++) {
	            if (i == csrRows.size() - 1) {
	                for (j = csrRows[i]-csrRows[0]; j < csrData.size(); j++) {
	                    //std::cout << "last row, " << i << ", " << j << ", rowSize = " << csrRows.size() << ", start = "
	                    //          << csrRows[i] << ", end = " << csrData.size() << ", csrData[" << j << "] = " << csrData[j]
	                    //          << std::endl;
	                    result[i] += csrData[j] * (double)denseVec[i];
	                }
	            } else {
	                for (j = csrRows[i]-csrRows[0]; j < csrRows[i + 1]-csrRows[0]; j++) {
	                    //std::cout << i << ", " << j << ", result size = " << result.size() << ", rowSize = " << csrRows.size()
	                    //          << " denseVec size = "
	                    //          << denseVec.size() << ", csrData size = " << csrData.size() << " , start = " << csrRows[i]
	                    //          << ", end = " << csrRows[i + 1] << ", csrData[" << j << "] = " << csrData[j] << std::endl;
	                    result[i] += csrData[j] * (double)denseVec[i];
	                }
	            }
	        }
	    }
	}
}

void csrSpMV::masterOnlySpMV(controlData control) {
    //convert sparse matrix from Matrix Market to Compressed Sparse Row format
    std::vector <Element> elements; //holds each matrix element as read from the Matrix Market format file
    int tempRow, tempCol;
    double tempData;

    // Read in sparse matrix saved in Matrix Market Format
    std::ifstream infile(control.matrixFile);
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
                        control.rowCount = std::stoi(token);
                    } else {
                        control.colCount = std::stoi(token);
                    }

                    i++;
                }

                control.nonZeros = std::stoi(line);
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
    printf("%d rows, %d cols, and %d non-zeros\n", control.rowCount, control.colCount, control.nonZeros);

    // add to csr vectors for use as CSR Format
    csrCols.resize(control.nonZeros);
    csrData.resize(control.nonZeros);

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

    if (control.rowCount != csrRows.size()){
        std::cout << "Actual row count does NOT match reported rowCount" << std::endl;
    }


    //std::cout << "Creating Dense Vector" << std::endl;
    if (!control.myId) {
        //std::cout << "Populating Dense Vector" << std::endl;
        for (int i = 0; i < control.rowCount; i++) {
            denseVec.push_back(1.0);
        }
    }

    std::cout << "Performing Master Only SpMV" << std::endl;
    result.resize(control.rowCount);

    std::cout << "rowcount = " << control.rowCount << ", csrRows.size() = " << csrRows.size() << std::endl;
    for (int i = 0; i < control.rowCount; i++) {
        //std::cout << "i = " << i << std::endl;
        double temp = 0.0;
        if (i != control.rowCount - 1) {
            if (control.colMajor) { //col major order selected
                for (int j = csrRows[i]; j < csrRows[i + 1]; j++) {   // go to end of current row
                    // entire row is multiplied by a single dense vector element
                    //std::cout << "i = " << i << " j = " << j << ", origin_row["<< i << "] = " << origin_row[i] << std::endl;
                    temp += csrData[j] * (double)denseVec[i];
                }
            } else {    // row major order selected
                for (int j = csrRows[i]; j < csrRows[i + 1]; j++) {   // go to end of current row
                    //std::cout << "i = " << i << " j = " << j << ", origin_row["<< i << "] = " << origin_row[i] << std::endl;
                    temp += csrData[j] * (double)denseVec[csrCols[j]];
                }
            }
        } else {
            if (control.colMajor) { //col major order selected
                for (int j = csrRows[i]; j < csrData.size(); j++) {  // go to end of data vector
                    // entire row is multiplied by a single dense vector element
                    //std::cout << "i = " << i << " j = " << j << ", origin_row["<< i << "] = " << origin_row[i] << std::endl;
                    temp += csrData[j] * (double)denseVec[i];
                }
            } else {    // row major order selected
                for (int j = csrRows[i]; j < csrData.size(); j++) {  // go to end of data vector
                    //std::cout << "i = " << i << " j = " << j << ", origin_row["<< i << "] = " << origin_row[i] << std::endl;
                    temp += csrData[j] * (double)denseVec[csrCols[j]];
                }
            }
        }
        result[i] = temp;
    }
}

