INC="-I${MKLROOT}/include -I${IZHROOT}/include"
LIB="-L${MKLROOT}/lib/intel64 -L${IZHROOT}/lib"
CFLAGS="gcc"
OPT="-g -O3 -Wall -std=c11"
LDFLAGS="-lizh_neuralnet -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm"

echo "compile run simulation.c"
$CFLAGS $OPT -o run_simulation.out run_simulation.c $INC $LIB $LDFLAGS