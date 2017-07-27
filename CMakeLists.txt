cmake_minimum_required(VERSION 3.5)
project(Distrubuted_SpMV)

set(CMAKE_CXX_STANDARD 11)
SET( CMAKE_C_COMPILER "/usr/bin/mpicc" )
SET( CMAKE_CXX_COMPILER "/usr/bin/mpic++" )

set(SOURCE_FILES main.cpp)
add_executable(Distrubuted_SpMV ${SOURCE_FILES})