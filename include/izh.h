#ifndef IZH
#define IZH

#include "ntk.h"

#define _block_size;
extern double _dt;


typedef enum _SYN_TYPE {
    BACKGROUND = 0,
    NO_DELAY,
    DELAY,
} SYN_TYPE ;


typedef struct _neuron_t
{
    int num_cells;
    double *v, *u, *ic;
    // izhikevich neuron model parameters
    double *a, *b, *c, *d; 
    // spike info
    int *num_spk;   
    int **t_fired; // fired time step for each neurons
    int *id_fired_neuron;

} neuron_t;


typedef struct _syn_t
{
    int num_pres;
    int num_syns;
    SYN_TYPE type;
    
    // size = num_syns
    int *id_pre_neuron;
    int *id_post_neuron;

    double *r, **ptr_r;
    double *veq, *weight;
    double *inv_tau;

    double **ptr_vpost;
    double **ptr_ipost;
    
    // delay parameter
    int *delay;
    
    // for background synapse
    double *p_fire;

} syn_t;

void init_cell_vars(neuron_t *cells, int num_cells, double cell_params[][4], int *cell_types);
void init_syn_vars(syn_t *syns, int num_pres, SYN_TYPE type, ntk_t *ntk, double cell_veq[], double cell_tau[], double *vpost, double *ipost);

void add_isyn_bck(syn_t *syns);
void add_isyn(syn_t *syns);
void add_isyn_delay(syn_t *syns);

void update_neurons(neuron_t *cells, int nstep, int *id_fired_neuron);
void update_syns_no_delay(syn_t *syns, int *id_fired_neuron);
void update_syns_delay(syn_t *syns); // need to udpate

int *gen_bck_spike(syn_t *bck_syns);

double *f_dv(double *v, void *arg_neuron, void *arg_null);
double *f_du(double *u, void *arg_neuron, void *arg_null);
double *f_dr_syns_no_delay(double *r, void *arg_syn, void *arg_fired);
double *f_dr_syns_delay(double *r, void *arg_syn, void *arg_fired);
double *solve_deq_using_rk4(double* (*f) (double*, void*, void*), int N, double *x, void *arg1, void *arg2);
void append_spike(int nstep, int *num_spk, int **t_spk);
double *read_ptr(int num_x, double **ptr_x);
void get_avg(int num_x, int **ptr_x);

void free_neurons(neuron_t *cells);
void free_syns(syn_t *syns);

#endif