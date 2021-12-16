#!/bin/bash

# usage: compile {fname} {output f name (optional)}
if [ -z "$2" ]; then
    val="./a.out"
else
    val=$2
fi

# gcc -o lib/Izh.o -c lib/Izh.c -O2 -Wall -I/opt/OpenBLAS/include -pthread -L/opt/OpenBLAS/lib -lm -lpthread -lgfortran -lopenblas
# gcc -o lib/Izh.o -c lib/Izh.c -fopenmp -O2 -Wall -I${MKLROOT}/include ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -lpthread -lm -ldl
gcc -o lib/Izh.o -g -I${MKLROOT}/include -c lib/Izh.c -lpthread -lm -ldl -fopenmp # ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -lpthread -lm -ldl -fopenmp
gcc -o lib/parson.o -c lib/parson.c
gcc -o lib/mt64.o -c lib/mt64.c
gcc -o lib/ntk.o -c lib/ntk.c
gcc -o tmp.o -c $1 -O2 -Wall -lm

# OPENBLAS
# gcc -Wall -O2 ./tmp.o lib/Izh.o lib/parson.o lib/mt64.o -o $2 -L/opt/OpenBLAS/lib -lm -lpthread -lgfortran -lopenblas
# MKL
gcc -g -Wall -O2 ./tmp.o lib/Izh.o lib/parson.o lib/mt64.o lib/ntk.o -o $2 -I${MKLROOT}/include ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -lpthread -lm -ldl -fopenmp

# gcc -Wall -O2 ./test/test.o ./lib/lib_Izh.o ./lib/mt64.o -o $2 -L/opt/OpenBLAS/lib -lm -lpthread -lgfortran -lopenblas

# gcc -o ./test/test.o -c $1 -O2 -Wall -lm
# gcc -Wall -O2 ./test/test.o ./lib/lib_Izh.o ./lib/mt64.o -o $2 -L/opt/OpenBLAS/lib -lm -lpthread -lgfortran -lopenblas

# gcc -I/opt/OpenBLAS/include -pthread -O2 -Wall $1 -o $2 -L/opt/OpenBLAS/lib -lm -lpthread -lgfortran -lopenblas
