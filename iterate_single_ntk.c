#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "parson.h"
#include "Izh.h"
#include "ntk.h"
#include "mt64.h"
#include "utils.h"

extern double gRatio; // from Izh.c
extern double CELL_TAU[_n_types];


typedef struct _simulinfo_t {

    int num_cells;
    double cell_ratio[_n_types];
    int *cell_types;

    double g_syn[_n_types][_n_types];
    double p_syn[_n_types][_n_types];

    int num_bck;
    double p_fire_bck, tau_bck;
    double g_bck[_n_types]; // coupling strength for each cell types
    double p_bck[_n_types]; // connection probability
    
    // name tag of the output file
    char tag[100];

} simulinfo_t;


int main(int argc, char **argv)
{

    init_genrand64(time(NULL));

    simulinfo_t info;

    info.num_cells = 100;
    
    info.cell_ratio[0] = 0.8;
    info.cell_ratio[1] = 0.2;
    reset2d(info.g_syn, 0); // TODO: change 
    reset2d(info.p_syn, 0.1);
    for (int i=0; i<_n_types; i++) {info.p_syn[1][i] = 0.4;} // 필요할까?

    info.num_bck   = 1000;
    info.p_fire_bck = 10 * _dt;
    info.tau_bck   = 5;
    reset1d(info.g_bck, 0);
    reset1d(info.p_bck, 0.1);

    info.g_bck[0] = 0.126; // check parameter
    info.g_bck[1] = 0.1;   // check parameter

    info.p_bck[0] = 0.1;
    info.p_bck[1] = 0.1;


}


void reset2d(double arg[_n_types][_n_types], double x)
{
    for (int i=0; i<_n_types; i++){
        for (int j=0; j<_n_types; j++){
            arg[i][j] = x;
        }
    }
}


void reset1d(double arg[_n_types], double x)
{
    for (int i=0; i<_n_types; i++){
        arg[i] = x;
    }
}


// void open_file_with_tag(FILE **fid, char tag[], char typename[])
// {
//     char fname[100];
//     sprintf(fname, "%s_%s.txt", tag, typename);
//     *fid = fopen(fname, "w");
// }
