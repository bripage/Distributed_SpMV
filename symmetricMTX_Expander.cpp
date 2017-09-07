//
// Created by brianpage on 9/7/17.
//
#include <iostream>
#include <string>
#include <istream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* strtod */

bool sortcol(const std::vector<double>& v1, const std::vector<double>& v2 ) {
    if (v1[1] == v2[1]){
        if (v1[0] < v2[0]) {
            return true;
        } else {
            return false;
        }
    } else if ( v1[1] < v2[1]){
        return true;
    } else {
        return false;
    }
    //return v1[1] < v2[1];
}

int main(int argc, char *argv[]) {
    std::string argTemp;
    std::string inputDataFile, outputDataFile;

    for (int i = 1; i < argc; i= i+2){
        argTemp = argv[i];
        if (argTemp == "-i"){
            // load data file.
            inputDataFile = argv[i+1];
        } else if (argTemp == "-o") {
            // save converted data file.
            outputDataFile = argv[i + 1];
        } else if (argTemp == "--help"){
            std::cout << "Usage: Symmetric MTX Expander [OPTION] <argument> ..." << std::endl << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << " -i <file>";
            std::cout << " -o <file>";
            exit(0);
        } else {
            printf("%s Is not a valid parameter. EXITING!\n", argv[i]);
            exit(0);
        }
    }

    std::vector<std::vector <double> > expandedMatrix;
    std::vector<double> tempNNZ;

    // Read in sparse matrix saved in Matrix Market Format
    std::ifstream infile(inputDataFile);
    if (!infile) {
        std::cout << "FAILED TO OPEN FILE!" << std::endl;
        exit(1);
    }

    int tempRow, tempCol, rowCount, colCount, nonZeros;
    double tempData;
    std::string matrixPropteries;
    std::string line;

    int i = 0, j = 0, previousRow = -1;
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
                        rowCount = std::stof(token);
                    } else {
                        colCount = std::stof(token);
                    }

                    i++;
                }

                nonZeros = std::stof(line);
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
                //std::cout << std::setprecision(20) << tempData << std::endl;

                tempNNZ.push_back(tempRow);
                tempNNZ.push_back(tempCol);
                tempNNZ.push_back(tempData);
                //std::cout << tempNNZ[2] << std::endl;
                expandedMatrix.push_back(tempNNZ);
                tempNNZ.clear();

                if (tempRow!= tempCol){
                    tempNNZ.push_back(tempCol);
                    tempNNZ.push_back(tempRow);
                    tempNNZ.push_back(tempData);
                    //std::cout << tempNNZ[2] << std::endl;
                    expandedMatrix.push_back(tempNNZ);
                    tempNNZ.clear();
                }
            }
        } else {
            previousLineCommented = true;
        }
    }


    // sort non zeros by column
    sort(expandedMatrix.begin(), expandedMatrix.end(), sortcol);

    std::ofstream outFile;
    outFile.open(outputDataFile, std::ofstream::out); // open the output file for writing the converted data set

    std::string temp;
    if (outFile.is_open()){
        std::cout << "Output file opened successfully" << std::endl;

        matrixPropteries = "";
        matrixPropteries += std::to_string(rowCount);
        matrixPropteries += " ";
        matrixPropteries += std::to_string(colCount);
        matrixPropteries += " ";
        matrixPropteries += std::to_string(expandedMatrix.size());
        matrixPropteries += "\n";

        outFile << matrixPropteries;

        // for each non zero, output it to file
        //outFile << std::fixed << std::setprecision(20);
        //outFile.flags (std::ios::scientific);
        //outFile.precision (std::numeric_limits<double>::digits10 + 1);
        for (int k = 0; k < expandedMatrix.size(); k++){
            outFile << expandedMatrix[k][0];
            outFile << " ";
            outFile << expandedMatrix[k][1];
            outFile << " ";
            outFile << std::setprecision(std::numeric_limits<long double>::digits10 -1) << expandedMatrix[k][2];
            outFile << std::endl;
            //std::cout << temp << std::endl;

            //outFile << std::setprecision(std::numeric_limits<long double>::digits10) << temp;
        }
    } else {
        std::cout << "Error opening output file";
        exit(0);
    }
    std::cout << "closing " << outputDataFile << std::endl;
    outFile.close(); // close the output file

    return 0;
}

