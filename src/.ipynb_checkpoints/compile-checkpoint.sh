#!/bin/bash

# usage: compile {fname} {output f name (optional)}
if [ -z "$2" ]; then
    val="./a.out"
else
    val=$2
fi

gcc -o ../lib/lib_Izh.o -c ../lib/lib_Izh.c -O2 -Wall -I/opt/OpenBLAS/include -pthread -L/opt/OpenBLAS/lib -lm -lpthread -lgfortran -lopenblas
gcc -o ../lib/parson.o -c ../lib/parson.c
gcc -o ../lib/mt64.o -c ../lib/mt64.c
gcc -o ./tmp.o -c $1 -O2 -Wall -lm

gcc -Wall -O2 ./tmp.o ../lib/lib_Izh.o ../lib/parson.o ../lib/mt64.o -o $2 -L/opt/OpenBLAS/lib -lm -lpthread -lgfortran -lopenblas

# gcc -Wall -O2 ./test/test.o ./lib/lib_Izh.o ./lib/mt64.o -o $2 -L/opt/OpenBLAS/lib -lm -lpthread -lgfortran -lopenblas

# gcc -o ./test/test.o -c $1 -O2 -Wall -lm
# gcc -Wall -O2 ./test/test.o ./lib/lib_Izh.o ./lib/mt64.o -o $2 -L/opt/OpenBLAS/lib -lm -lpthread -lgfortran -lopenblas

# gcc -I/opt/OpenBLAS/include -pthread -O2 -Wall $1 -o $2 -L/opt/OpenBLAS/lib -lm -lpthread -lgfortran -lopenblas
