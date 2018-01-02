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

void csrSpMV::rebase(int colAdjustment){
	int firstRowPtr = csrRows[0];
	for (int i = 0; i < csrRows.size(); i++){
		csrRows[i] -= firstRowPtr;
	}

	for (int i = 0; i < csrCols.size(); i++){
		//std::cout << "csrCols[" << i << "] = " << csrCols[i] << " - " << colAdjustment << std::endl;
		csrCols[i] -= colAdjustment;
	}
}

void csrSpMV::nodeSpMV(controlData control, std::vector <double>& result) {

    //std::this_thread::sleep_for(std::chrono::milliseconds(control.myId * 2000));
    //std::cout << "myID = " << control.myId << ", data length = " << csrData.size() << std::endl;

    if (csrData.size() > 0) {
		int ompThreadId, start, end, i, j;

		#pragma omp parallel num_threads(control.ompThreads) shared(result) private(ompThreadId, start, end, i, j)
	    {
		    ompThreadId = omp_get_thread_num();

	        for (i = 0; i < csrRows.size(); i++) {
	            if (i == csrRows.size() - 1) {
	                for (j = csrRows[i]-csrRows[0]; j < csrData.size(); j++) {
	                    result[i] += csrData[j] * (double)denseVec[i];
	                }
	            } else {
	                for (j = csrRows[i]-csrRows[0]; j < csrRows[i + 1]-csrRows[0]; j++) {
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

                if (!(tempElement.row == 0 && tempElement.col == 0 && tempElement.data == 0.0)) {
                    elements.push_back(tempElement);
	                //std::cout << "tempdata = " << tempElement.data << std::endl;
                }
            }
        } else {
            previousLineCommented = true;
        }
    }
    printf("%d rows, %d cols, and %d non-zeros\n", control.rowCount, control.colCount, control.nonZeros);
/*
	std::cout << "Rows before sorting: " << std::endl;
	for (int i = 0; i < elements.size(); i++) {
		std::cout << elements[i].row << ",";
	}
	std::cout << std::endl;
*/
	// sort the vector of elements by row, and then each row, based on column
	std::stable_sort(elements.begin(), elements.end(), sortByCol);
/*
	std::cout << "Rows after sorting: " << std::endl;
	for (int i = 0; i < elements.size(); i++) {
		std::cout << elements[i].row << ",";
	}
	std::cout << std::endl;
*/
    // add to csr vectors for use as CSR Format
	int previousRow = -1;
    for (int k = 0; k < elements.size(); k++) {
	    //make sure all columns have a row pointer for the current elements row. If not create it!
	    if (elements[k].row == previousRow) {
		    // Do nothing, row is already present in list
	    } else {
		    // Add new row pointer, as it has not been seen previously
		    csrRows.push_back(csrData.size());
	    }
        csrCols.push_back(elements[k].col);
        csrData.push_back(elements[k].data);
	    //std::cout << "csrData[" << k << "] = " << csrData[k] << std::endl;

	    previousRow = elements[k].row;
    }


    if (control.rowCount != csrRows.size()){
        std::cout << "Actual row count does NOT match reported rowCount" << std::endl;
    }

/*
	std::cout << "Rows: " << std::endl;
	for (int i = 0; i < control.rowCount; i++) {
		std::cout << csrRows[i] << ",";
	}
	std::cout << std::endl;
	std::cout << std::endl << "Master Only Distribution of NonZero Elements" << std::endl;
	for (int j = 0; j < csrData.size(); j++) {
		std::cout << csrData[j] << ",";
	}
	std::cout << std::endl;
*/


    // Fill the dense vector with preliminary data for use in the Master Only SpMV calculation
    if (!control.myId) {
        for (int i = 0; i < control.rowCount; i++) {
            denseVec.push_back(1.0);
        }
    }

    std::cout << "Performing Master Only SpMV" << std::endl;
    result.resize(control.rowCount);

    std::cout << "rowcount = " << control.rowCount << ", csrRows.size() = " << csrRows.size() << std::endl;


	for (int i = 0; i < control.rowCount; i++) {
		double temp = 0.0;
		if (i != control.rowCount - 1) {
			for (int j = csrRows[i]; j < csrRows[i + 1]; j++) {   // go to end of current row
				// entire row is multiplied by a single dense vector element
				std::cout << i << "," << j << " - " << csrData[j] << " * " << denseVec[csrCols[j]] << std::endl;
				temp += csrData[j] * (double) denseVec[csrCols[j]];
			}
		} else {
			for (int j = csrRows[i]; j < csrData.size(); j++) {  // go to end of data vector
				// entire row is multiplied by a single dense vector element
				std::cout << i << "," << j << " - " << csrData[j] << " * " << denseVec[csrCols[j]] << std::endl;
				temp += csrData[j] * (double) denseVec[csrCols[j]];
			}
		}
		result[i] = temp;
	}
}

