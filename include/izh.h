#ifndef _IZH
#define _IZH

#include "ntk.h"

// #define _dt 0.005
#define _block_size 500
double _dt;

typedef struct _neuron_t
{
    int num_cells;
    // dynamic variables
    double *v, *u, *ic;
    // parameters
    double *a, *b, *c, *d;
    double mu_wn, std_wn; // white noise

    // spike info
    int *num_spk; // index to append the new spike
    int **t_spk;  // activated steps for each neuron
    int *id_fire; // activated cell id

} neuron_t;


typedef struct _syn_t
{   
    // assume that network is sparse
    // indicators
    int is_delay_on; // indicates the synapse exists the delay

    int num_cells; // the # of the nodes in the network
    int num_syns;  // the # of the synapses
    int num_pre;

    int *id_pre_neuron; // presynaptic neuron id
    int *id_post_neuron; // postsynaptic neuron id
    
    double R;
    double *r, *weight;
    double *veq, *inv_tau;
    double **ptr_vpost; // pointer to the voltage of post-neuron
    double **ptr_ipost;

    /*** delay_on == 1 ***/
    int *delay;
    int *id_pre; // (num_syns, )
    int *id_t_spk; // (num_syn, )

    /*** delay_on == 0 ***/
    int *n_pre2syn; // (num_pre, )
    int **id_pre2syn; // (num_pre, num_edges)

} syn_t;


typedef struct _bcksyn_t
{   
    // use matrix-multiplication to update synaptic state
    // (because bck_syn is very densed network)

    int num_cells;
    int num_bck;
    int num_pre;

    // assume that background synapse is always the excitatory
    double R;
    double *r, *weight; // size = num_cells
    double *inv_tau; 
    double *fspk;
    double *syn_act; // if the background synapse activated, 1

} bcksyn_t;


typedef struct _writer_t
{
    char tag[50];
    int N;
    FILE *fv, *fu, *fr, *fic, *fid;

} writer_t;


// init functions
void init_random_stream(long int seed);
void init_cell_vars(neuron_t *cells, int num_cells, int *cell_types);
void init_syn_vars(syn_t *syns, int num_cells, int is_delay_on, double *delay, double *vpost, double *ipost, ntk_t *ntk);
void init_bcksyn_vars(bcksyn_t *bck_syns, int num_cells, int num_bck, double tau_bck, ntk_t *ntk);
void map_pre2syn(syn_t *syns, ntk_t *ntk);

// update functions
void update_no_delay(int nstep, neuron_t *cells, syn_t *syns, bcksyn_t *bck_syns);
void update_delay(int nstep, neuron_t *cells, syn_t *syns, bcksyn_t *bck_syns);
void update_cells(int nstep, neuron_t *cells);
void update_syns_no_delay(int nstep, syn_t *syns, void *id_fire_tmp);
void update_syns_delay(int nstep, syn_t *syns, int **t_spk);
void update_bcksyns(bcksyn_t *bck_syns);
void get_bck_input(syn_t *bck_syns);
void add_syn_current(syn_t *syns, neuron_t *cells);
void add_bcksyn_current(bcksyn_t *bck_syns, neuron_t *cells);
void append_spike(int nstep, int *id, int **t_spk);

// related to solve differential equation
void solve_deq_using_rk4(void (*f) (double*, double*, void*, void*), int N, double *dx, double *x, void *arg1, void *arg2);
void f_dv(double *dv, double *v, void *arg_cells, void *arg_null);
void f_du(double *du, double *u, void *arg_cells, void *arg_null);
void f_dr_delay(double *dr, double *r, void *arg_syns, void *arg_syn_act);
void f_dr_no_delay(double *dr, double *r, void *arg_syns, void *arg_id_fire);
void f_dr_bck(double *dr, double *r, void *arg_syns, void *arg_syn_act);

// mathematical functions
double get_avg(int num_x, double *x);

// Kuramoto order parameter
void get_Kuramoto_order_params(int len, neuron_t *cells, double *rK, double *psiK);
void get_spike_phase(int n_spk, int nmax, int *nsteps, double *phase);

// utilities
void free_cells(neuron_t *cells);
void free_syns(syn_t *syns);
void free_bcksyns(bcksyn_t *bck_syns);
void reset_spike(neuron_t *cells);

// save data
void save_env(char fname[100], int num_cells, int *cell_types, int num_syns, int *id_presyns, double tmax);
void write_dat(FILE *fid, int N, double *x);
void write_cell_fire(FILE *fid, int *id_fire);

#endif