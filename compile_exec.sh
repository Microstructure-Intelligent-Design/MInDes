#! /bin/bash
rm -f *.o
echo "> Program is compiling... Please wait..."
g++ -std=c++11 -fopenmp -static -c -w solvers/Solver_Allen_Cahn.cpp -o Solver_Allen_Cahn.o
g++ -std=c++11 -fopenmp -static -c -w solvers/Solver_Cahn_Hilliard.cpp -o Solver_Cahn_Hilliard.o
g++ -std=c++11 -fopenmp -static -c -w solvers/Solver_custom.cpp -o Solver_custom.o
g++ -std=c++11 -fopenmp -static -c -w solvers/Solver_Fluid.cpp -o Solver_Fluid.o
g++ -std=c++11 -fopenmp -static -c -w solvers/Solver_Mechanics.cpp -o Solver_Mechanics.o
g++ -std=c++11 -fopenmp -static -c -w solvers/Solver_poisson_equation.cpp -o Solver_poisson_equation.o
g++ -std=c++11 -fopenmp -static -c -w solvers/Solver_Temperature.cpp -o Solver_Temperature.o
g++ -std=c++11 -fopenmp -static -c -w solvers/Solvers.cpp -o Solvers.o
g++ -std=c++11 -fopenmp -static -c -w MInDes_EXE.cpp -o MInDes_EXE.o
g++ -std=c++11 -fopenmp -static -o MInDes *.o -Llib/linux/lib -lfftw3
rm -f *.o
echo "> MInDes has been builded"
