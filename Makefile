CC = gcc
CFLAGS = -Wall -O2 -std=c11 -g 
# LDFLAGS = -lm -fopenmp -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
# LIB = -L${MKLROOT}/lib/intel64 
LDFLAGS = -lizh_neuralnet -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
LIB = -L${MKLROOT}/lib/intel64 -L${IZHROOT}/lib
INC = -I${MKLROOT}/include -I${IZHROOT}/include


TARGET = run_single.out

# need to update to LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/mkl/lib/intel64

.PHONY: main clean

main:
	$(CC) $(CFLAGS) $(TARGET:.out=.c) -o $(TARGET) $(INC) $(LIB) $(LDFLAGS)

clean:
	rm -f $(TARGET)