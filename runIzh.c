#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "lib/Izh.h"
#include "lib/ntk.h"
#include "lib/mt64.h"
#include "lib/easyIzh.h"
#include "omp.h"



int main()
{
    init_genrand64(100);

    simul_info_t info;
    empty_simulvar(&info);
    /*** SET PARAMETERS ***/
    info.num_cells = 100;

    // cell types
    info.cell_ratio[0] = 0.8;
    info.cell_ratio[1] = 0.2;

// synapse 
    double g, p;

    info.g_syn[0][0] = 3; // e->e
    info.g_syn[0][1] = 3; // e->i
    info.g_syn[1][0] = 16; // i->e
    info.g_syn[1][1] = 16; // i->i

    info.p_syn[0][0] = 0.1;
    info.p_syn[0][1] = 0.1;
info.p_syn[1][0] = 0.1;
    info.p_syn[0][1] = 0.1;

    // background input
    info.num_bck = 2000;
    info.p_fire_bck = 3 * _dt; // kHz
    info.tau_bck = 5;

    g = 0.6;
    info.g_bck[0] = g;
    info.g_bck[1] = g;

    p = 0.01;
    info.p_bck[0] = p;
    info.p_bck[1] = p;

    // time
    info.tmax = 2000;
    info.amp_sq[0] = 5;
    info.amp_sq[1] = 1;
    
    // run simulation
    run_simulation(&info);

}
