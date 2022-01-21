INC="-I${MKLROOT}/include -I${IZHROOT}/include"
LIB="-L${MKLROOT}/lib/intel64 -L${IZHROOT}/lib"
CFLAGS="gcc"
OPT="-g -O3 -Wall -std=c11"
LDFLAGS="-lizh_neuralnet -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm"

echo "$1 -> $2"
$CFLAGS $OPT -o $2 $1 subfunction.o $INC $LIB $LDFLAGS