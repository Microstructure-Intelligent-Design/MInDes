#! /bin/bash
rm -f *.o
echo "> Program is compiling... Please wait..."
g++ -std=c++11 -fopenmp -static -c -w solvers/*.cpp
g++ -std=c++11 -fopenmp -static -c -w MInDes_EXE.cpp
g++ -std=c++11 -fopenmp -static -o MInDes *.o -Llib/linux/lib -lfftw3
rm -f *.o
echo "> MInDes has been builded"
