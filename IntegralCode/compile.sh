# Compile Special Mathfunction library, create library
gcc  -std=c99 -Ofast -fassociative-math -flto -march=native -funroll-loops -fPIC  -c SpecialMath.c  -o spmath.o

# Compile Cython functions
#g++ -Ofast -flto -march=native -fopenmp -funroll-loops -fPIC -shared  -lrt -c iCodeCore.cpp
gcc -std=c99 -Ofast -flto -march=native -fopenmp -funroll-loops -fPIC -shared  -lrt -c iCodeCore.c

# Create cython file, compile to python module
#cython --cplus iCodeModule.pyx
cython iCodeModule.pyx

gcc -shared -fPIC -O3 -Wall -fno-strict-aliasing  -I/usr/include/python2.7 -o iCodeModule.so iCodeModule.c  iCodeCore.o spmath.o -lgomp -lm 

# Clean-up temporary files
rm *.o
#rm iCodeModule.c
rm -rf build

