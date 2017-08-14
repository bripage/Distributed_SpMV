#include <iostream>
#include <string>
#include <istream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <numeric>

int main(int argc, char *argv[]) {
    std::string argTemp;
    std::string inputDataFile, outputDataFile;
    int floatingPointOps = 0;
    int totalMPITasks = 0;

    for (int i = 1; i < argc; i= i+2){
        argTemp = argv[i];
        if (argTemp == "-i"){
            // load data file.
            inputDataFile = argv[i+1];
        } else if (argTemp == "-o") {
            // save converted data file.
            outputDataFile = argv[i + 1];
        } else if (argTemp == "--fops"){
            floatingPointOps = atoi(argv[i + 1]);
        } else if (argTemp == "--tasks"){
            totalMPITasks = atoi(argv[i + 1]);
        } else if (argTemp == "--help"){
            std::cout << "Usage: SpMVDataConvert [OPTION] <argument> ..." << std::endl << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << " -i <file>" << std::setw(15) << "" << R"(Load data file <file>)" << std::endl;
            std::cout << " -o <file>" << std::setw(15) << "" << R"(Save converted data file as <file>)" << std::endl;
            exit(0);
        } else {
            printf("%s Is not a valid parameter. EXITING!\n", argv[i]);
            exit(0);
        }
    }


    std::vector<double> distributionTimes, SmPVTimes, totalTimes;
    std::vector<std::vector <double> > combinedData;

    // Read in sparse matrix saved in Matrix Market Format
    std::ifstream infile(inputDataFile);
    if (!infile) {
        std::cout << "FAILED TO OPEN FILE!" << std::endl;
        exit(1);
    }
    std::string line;

    int i = 0, j = 0;
    int processCount = -1, threadCount = -1;
    double distSum = 0.0, spmvSum = 0.0, totalSum = 0.0;

    std::string headerLine = "MPI Processes,Nodes,Processes per Node,OpenMP Threads per process,Floating Point Ops," \
            "Distribution Time,SpMV Time,Total Time, ,Best Distribution Time,Best SpMV Time,Best Total Time,Best GFlops," \
            ",Worst Distribution Time,Worst SpMV Time,Worst Total Time,Worst GFlops, ,Avg Distribution Time," \
            "Avg SpMV Time,Avg Total Time,Avg GFlops\n";
    std::string blankLine= " , , , , , , , , , , , , \n";
    std::string dataLine;
    std::ofstream outFile;
    outFile.open(outputDataFile, std::ofstream::out); // open the output file for writing the converted data set

    if (outFile.is_open()){
        std::cout << "Output file opened successfully" << std::endl;

        while (std::getline(infile, line)) {
            if (line.find("task") != std::string::npos) {
                if (processCount == -1 || threadCount == -1) {
                    size_t pos = 0;
                    std::string token;

                    if (line.find("node") != std::string::npos) {
                        pos = line.find(' ');
                        token = line.substr(0, pos);

                        processCount = atoi(token.c_str());
                    } else if (line.find("threads") != std::string::npos) {
                        pos = line.find(' ');
                        token = line.substr(0, pos);

                        threadCount = atoi(token.c_str());
                    }
                } else {
                    if (!distributionTimes.empty() && !SmPVTimes.empty() && !totalTimes.empty()) {
                        if (distributionTimes[0] == 0 && SmPVTimes[0] == 0 && totalTimes[0] == 0){
                            distributionTimes.erase(distributionTimes.begin());
                            SmPVTimes.erase(SmPVTimes.begin());
                            totalTimes.erase(totalTimes.begin());
                        }
                        std::sort(distributionTimes.begin(), distributionTimes.end());
                        std::sort(SmPVTimes.begin(), SmPVTimes.end());
                        std::sort(totalTimes.begin(), totalTimes.end());

                        outFile << headerLine;
                        //outFile << std::endl;

                        std::vector<double> temp;

                        for (int i = 0; i < distributionTimes.size(); i++) {
                            if (i == 0) {
                                dataLine = "";
                                dataLine += std::to_string(totalMPITasks);
                                dataLine += ", ,";
                                dataLine += std::to_string(processCount);
                                dataLine += ",";
                                dataLine += std::to_string(threadCount);
                                dataLine += ",";
                                dataLine += std::to_string(floatingPointOps);
                                dataLine += ",";
                                dataLine += std::to_string(distributionTimes[0]);
                                dataLine += ",";
                                dataLine += std::to_string(SmPVTimes[0]);
                                dataLine += ",";
                                dataLine += std::to_string(totalTimes[0]);
                                dataLine += ", ,";
                                dataLine += std::to_string(distributionTimes[0]);
                                dataLine += ",";
                                dataLine += std::to_string(SmPVTimes[0]);
                                dataLine += ",";
                                dataLine += std::to_string(totalTimes[0]);
                                dataLine += ",";
                                dataLine += std::to_string(((floatingPointOps * 2) / SmPVTimes[0]) / 1000000000);
                                dataLine += ", ,";
                                dataLine += std::to_string(distributionTimes[distributionTimes.size() - 1]);
                                dataLine += ",";
                                dataLine += std::to_string(SmPVTimes[SmPVTimes.size() - 1]);
                                dataLine += ",";
                                dataLine += std::to_string(totalTimes[totalTimes.size() - 1]);
                                dataLine += ",";
                                dataLine += std::to_string(
                                        ((floatingPointOps * 2) / SmPVTimes[SmPVTimes.size() - 1]) / 1000000000);
                                dataLine += ", ,";
                                dataLine += std::to_string(distSum / double(distributionTimes.size()));
                                dataLine += ",";
                                dataLine += std::to_string(spmvSum / double(SmPVTimes.size()));
                                dataLine += ",";
                                dataLine += std::to_string(totalSum / double(totalTimes.size()));
                                dataLine += ",";
                                dataLine += std::to_string(
                                        (floatingPointOps * 2) / (spmvSum / double(SmPVTimes.size())) /
                                        1000000000);
                                dataLine += "\n";

                                temp.push_back(totalMPITasks);
                                temp.push_back(processCount);
                                temp.push_back(threadCount);
                                temp.push_back(distributionTimes[0]);
                                temp.push_back(SmPVTimes[0]);
                                temp.push_back(totalTimes[0]);
                                temp.push_back(((floatingPointOps * 2) / SmPVTimes[0]) / 1000000000);
                                temp.push_back(distributionTimes[distributionTimes.size() - 1]);
                                temp.push_back(SmPVTimes[SmPVTimes.size() - 1]);
                                temp.push_back(totalTimes[totalTimes.size() - 1]);
                                temp.push_back(((floatingPointOps * 2) / SmPVTimes[SmPVTimes.size() - 1]) / 1000000000);
                                temp.push_back(distSum / double(distributionTimes.size()));
                                temp.push_back(spmvSum / double(SmPVTimes.size()));
                                temp.push_back(totalSum / double(totalTimes.size()));
                                temp.push_back((floatingPointOps * 2) / (spmvSum / double(SmPVTimes.size())) /
                                1000000000);
                                temp.push_back((((floatingPointOps * 2) / SmPVTimes[0]) / 1000000000) -
                                                       (((floatingPointOps * 2) / SmPVTimes[SmPVTimes.size() - 1])
                                                        / 1000000000));
                                temp.push_back((((floatingPointOps * 2) / SmPVTimes[0]) / 1000000000) /
                                                       (((floatingPointOps * 2) / SmPVTimes[SmPVTimes.size() - 1])
                                                        / 1000000000));

                                combinedData.push_back(temp);

                            } else if (i != 0) {
                                dataLine = " , , , , ,";
                                dataLine += std::to_string(distributionTimes[i]);
                                dataLine += ",";
                                dataLine += std::to_string(SmPVTimes[i]);
                                dataLine += ",";
                                dataLine += std::to_string(totalTimes[i]);
                                dataLine += " , , , , , , , , , , , , , , \n";
                            }

                            outFile << dataLine;
                        }

                        outFile << blankLine;
                        outFile << blankLine;

                        distributionTimes.clear();
                        SmPVTimes.clear();
                        totalTimes.clear();
                        distSum = 0;
                        spmvSum = 0;
                        totalSum = 0;
                    }
                }

                size_t pos = 0;
                std::string token;

                if (line.find("node") != std::string::npos){
                    pos = line.find(' ');
                    token = line.substr(0, pos);

                    processCount = atoi(token.c_str());
                } else if (line.find("threads") != std::string::npos){
                    pos = line.find(' ');
                    token = line.substr(0, pos);

                    threadCount = atoi(token.c_str());
                }
            } else {
                size_t pos = 0;
                std::string token;

                pos = line.find_first_of(',');
                token = line.substr(0, pos);
                double dist = atof(token.c_str());
                line.erase(0, pos + 1);

                pos = line.find_first_of(',');
                token = line.substr(0, pos);
                double spmv = atof(token.c_str());
                line.erase(0, pos + 1);

                token = line.substr(0, pos);
                double total = atof(token.c_str());

                distributionTimes.push_back(dist);
                distSum += dist;
                SmPVTimes.push_back(spmv);
                spmvSum += spmv;
                totalTimes.push_back(total);
                totalSum += total;
            }
        }
        outFile << blankLine;
        outFile << blankLine;
        outFile << blankLine;

        outFile << "MPI Processes, Processes per Node, Threads per Process, ,Best Distribution Time, Best SpMV Time," \
                    "Best Total Time, Best GFlops, ,Worst Distribution Time, Worst SpMV Time, Worst Total Time," \
                    "Worst GFlops, ,Avg Distribution Time, Avg SpMV Time, Avg Total Time, Avg GFlops, ," \
                    "Performance Range, Performance Ratio (Best/Worst), \n";
        for (int i = 0; i < combinedData.size(); i++){
            std::string tempString = "";
            tempString += std::to_string(combinedData[i][0]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][1]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][2]);
            tempString += ", ,";
            tempString += std::to_string(combinedData[i][3]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][4]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][5]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][6]);
            tempString += ", ,";
            tempString += std::to_string(combinedData[i][7]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][8]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][9]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][10]);
            tempString += ", ,";
            tempString += std::to_string(combinedData[i][11]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][12]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][13]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][14]);
            tempString += ", ,";
            tempString += std::to_string(combinedData[i][15]);
            tempString += ",";
            tempString += std::to_string(combinedData[i][16]);
            tempString += ", \n";

            outFile << tempString;
        }
    } else {
        std::cout << "Error opening output file";
        exit(0);
    }
    std::cout << "closing " << outputDataFile << std::endl;
    outFile.close(); // close the output file

    return 0;
}
