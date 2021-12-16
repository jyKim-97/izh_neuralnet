#ifndef _EASY_IZH
#define _EASY_IZH

#include "Izh.h"

typedef struct _simul_info_t {

    int num_cells;
    double cell_ratio[_n_types];
    int *cell_types;

    double g_syn[_n_types][_n_types];
    double p_syn[_n_types][_n_types];

    int num_bck;
    double p_fire_bck, tau_bck;
    double g_bck[_n_types]; // coupling strength for each cell types
    double p_bck[_n_types]; // connection probability

    // option
    int is_delay_on;
    double *delay;

    double mu_noise, std_noise;
    double tmax;
    
    // square current input
    double amp_sq[_n_types];
    double t_sq_inp[2];

    // name tag of the output file
    char tag[100];

} simul_info_t;

void read_simulinfo(simul_info_t *info, char jsonfname[100]);
void run_simulation(simul_info_t *info);
void init_simulation(simul_info_t *info, neuron_t *cells, syn_t *syns, bcksyn_t *bck_syns);
void empty_simulvar(simul_info_t *info);
void gen_cell_type(int num_cells, double cell_ratio[], int *cell_types);
void square_current();

// for testing
void print_ntk_info(ntk_t *ntk, int n_print_node);
void print_post2syn_id(syn_t *syns, int n_print_node);
void print_variable(double *x, int n_print_node);

#endif