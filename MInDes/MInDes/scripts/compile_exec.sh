#! /bin/bash
rm -f *.o
echo "> Program is compiling... Please wait..."
g++ -std=c++17 -fopenmp -O3 -march=native -c -w src/modules/*.cpp &
g++ -std=c++17 -fopenmp -O3 -march=native -c -w src/modules/base/preprocess_modules/MicrostructureInit/*.c &
wait
g++ -std=c++17 -fopenmp -O3 -march=native -o MInDes *.o -Llib/fftw/linux/lib -lfftw3
rm -f *.o
echo "> MInDes has been builded"
