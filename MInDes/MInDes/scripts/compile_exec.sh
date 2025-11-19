#! /bin/bash
rm -f *.o
NPROC=$(nproc 2>/dev/null || echo 4)
echo "> Program is compiling... Please wait..."
g++ -std=c++17 -flto="$NPROC" -fopenmp -O3 -march=x86-64 -mtune=generic -c -w src/*.cpp &
g++ -std=c++17 -flto="$NPROC" -fopenmp -O3 -march=x86-64 -mtune=generic -c -w src/modules/*.cpp &
wait
g++ -std=c++17 -flto="$NPROC" -fopenmp -O3 -march=x86-64 -mtune=generic -o MInDes *.o -Llib/linux/lib -lfftw3 -lMInDeslib
rm -f *.o
echo "> MInDes has been built"
