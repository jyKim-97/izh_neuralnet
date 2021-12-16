CC = gcc
CFLAGS = -Wall -O2 -std=c11 -g 
LDFLAGS = -lm -fopenmp
INC = -I/opt/intel/mkl/include -Iinclude 
LIB_MKL = -L/opt/intel/mkl/lib/intel64/ -ldl -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread 

TARGET = runIzh.out
SINGLE_TARGET = run_single.out


# need to update to LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/mkl/lib/intel64

SRCS = include/Izh.c include/mt64.c include/ntk.c include/parson.c include/utils.c # lib/easyIzh.c 
OBJS = $(SRCS:.c=.o)
# OBJS_P = $(SRCS:.c=.o)

.PHONY: main test clean

main: $(OBJS)
	$(CC) $(CFLAGS) $(INC) -o include/main.o -c $(TARGET:.out=.c)
	$(CC) $(CFLAGS) $(INC) $(OBJS) include/main.o -o $(TARGET) $(LIB_MKL) $(LDFLAGS)


single: $(OBJS)
	$(CC) $(CFLAGS) $(INC) -o include/tmp.o -c $(SINGLE_TARGET:.out=.c)
	$(CC) $(CFLAGS) $(INC) $(OBJS) include/tmp.o -o $(SINGLE_TARGET) $(LIB_MKL) $(LDFLAGS)


# profiling: $(OBJS_P)
# 	$(CC) $(CFLAGS) -pg -o lib/tmp_p.o -c $(SINGLE_TARGET:.out=.c)
# 	$(CC) $(CFLAGS) -pg $(OBJS_P) lib/tmp_p.o -o $(SINGLE_TARGET:.out=_p.out) $(LIB_MKL) $(LDFLAGS)


$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) -o $@ -c $*$..c $(INC) $(LDFLAGS)



clean:
	rm -f lib/*.o
	rm -f $(TARGET)