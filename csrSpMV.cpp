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
	processData = objToCopy.processData;
	csrRows = objToCopy.csrRows;
	csrCols = objToCopy.csrCols;
	csrData = objToCopy.csrData;
}

csrSpMV::~csrSpMV() {
	processData.clear();
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
		csrCols[i] = csrCols[i] - colAdjustment;
	}
}

void csrSpMV::rebase_balanced(){
	int firstRowPtr = csrRows[0];
	for (int i = 0; i < csrRows.size(); i++){
		csrRows[i] -= firstRowPtr;
	}
}

void csrSpMV::nodeSpMV(controlData control, std::vector <double>& result) {
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

void csrSpMV::masterOnlySpMV(controlData control, std::vector<int>& seqDist) {
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
	                    result.resize(control.rowCount, 0.0);
	                    denseVec.resize(control.rowCount, 1.0);
                    } else {
                        control.colCount = std::stoi(token);
                    }

                    i++;
                }
                control.nonZeros = std::stoi(line);
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

	            if (!(tempData == 0.0 || tempData == 0)) {
	                result[tempRow] += tempData * denseVec[tempCol];
                    seqDist[tempRow]++;
                }
            }
        } else {
            previousLineCommented = true;
        }
    }
}

