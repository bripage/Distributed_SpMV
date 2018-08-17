main :
  mpic++ -std=c++11 -fopenmp -lm -O3 main.cpp csrSpMV.cpp distribution.cpp -o main

clean :
	rm  main
