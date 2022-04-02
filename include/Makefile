# compile static library
CC = gcc
CFLAGS = -Wall -O3 -std=c11 -g 
LDFLAGS = -lm -fopenmp -lmkl_sequential -lmkl_core -lpthread
INC = -I${MKLROOT}/include
LIB = -L${MKLROOT}/lib/intel64/

TARGET = ../lib/libizh_neuralnet.a

SRCS = izh.c mt64.c parson.c utils.c ntk.c writer.c
OBJS = $(SRCS:.c=.o)

.PHONY: clean

main: $(OBJS)
	ar -rc $(TARGET) -o $(OBJS)
#	ar -rcT libtest.a lizh_neuralnet


$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) -o $@ -c $*$..c $(INC) $(LIB) $(LDFLAGS)

clean:
	rm -f *.o
