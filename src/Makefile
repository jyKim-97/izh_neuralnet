CC = gcc
CFLAGS = -Wall -O2 -I/opt/OpenBLAS/include
LDFLAGS = -L/opt/OpenBLAS/lib -lm -lpthread -lgfortran -lopenblas
TARGET = ./run_Izh_random_input.c

SRCS = ../lib/mt64.c ../lib/lib_Izh.c ../lib/parson.c ./run_Izh_random_input.c
OBJS = ../lib/mt64.o ../lib/lib_Izh.o ../lib/parson.o ./run_Izh_random_input.o

$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) -o $@ -c $< $(LDFLAGS)

# main: $(OBJS)
# 	$(CC) $(CFLAGS) $(LDFLAGS) -o ./run.out