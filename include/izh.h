#ifndef IZH
#define IZH

#include "ntk.h"
#define MAX_TYPE 6

extern double _dt;




typedef enum _SYN_TYPE {
    BACKGROUND = 0x01 << 0,
    NO_DELAY = 0x01 << 1,
    DELAY = 0x01 << 2,
    STD = 0x01 << 3,
    STF = 0x01 << 4
} SYN_TYPE;


typedef enum _NET_TYPE {
    MEAN_DEG = 0,
    PROB
} NET_TYPE;


typedef struct _neuron_t {
    int num_cells;
    double *v, *u, *ic;
    // izhikevich neuron model parameters
    double *a, *b, *c, *d; 
    // spike info
    int *num_spk;  
    int **t_fired; // fired "time step" for each neurons
    int *id_fired;
    int nstep;
    int *types; // cell types

} neuron_t;


typedef struct _syn_t {
    int num_pres;
    int num_syns;
    SYN_TYPE type;
    
    // size = num_syns
    int *id_pre_neuron;
    int *id_post_neuron;

    double *r, *veq;
    double *inv_tau;
    double *weight;
    double **ptr_ipost, **ptr_vpost;

    // short-term plasticity method (depression)
    double *x, *z; // x+r+z=1
    double *u;
    double tau_r;
    double U, tau_facil;
    
    // delay parameter
    int *delay;
    int *id_exp; // explore id, use in delay
    
    // for background synapse
    double *p_fire;

} syn_t;


typedef struct _network_info_t {
    int num_cells;
    double cell_type_ratio[MAX_TYPE];
    double cell_params[MAX_TYPE][4]; // num_types * 4 or NULL

    int mean_degs[MAX_TYPE][MAX_TYPE];
    double psyns[MAX_TYPE][MAX_TYPE]; // pre (row) / post (col)
    double gsyns[MAX_TYPE][MAX_TYPE];
    double syn_tau[MAX_TYPE]; // num_types
    double syn_veq[MAX_TYPE]; // num_types

    int num_bck;
    double bck_type_ratio[MAX_TYPE];

    double pbck[MAX_TYPE][MAX_TYPE];
    double gbck[MAX_TYPE][MAX_TYPE];
    double bck_tau[MAX_TYPE];
    double bck_veq[MAX_TYPE];
    double frbck[MAX_TYPE]; // size = num_bck_types

    double tmax;
    double fs; // sampling rate
    double t_delay_m, t_delay_std;

    SYN_TYPE type;
    NET_TYPE type_ntk;

} network_info_t;


typedef struct _result_t {
    double fs; // sampling rate
    int num_times;
    double *t, *vm, *rk;
    int num_freqs;
    double *freq, *yf;
} res_t;


typedef struct _arg_t {
    double *x, *y;
    int len;
    int message;
} arg_t;


void init_random_stream(long int seed);
void init_cell_vars(neuron_t *cells, int num_cells, double cell_params[MAX_TYPE][4], int *cell_types);
void init_syn_vars(syn_t *syns, int num_pres, SYN_TYPE type, ntk_t *ntk, double cell_veq[], double cell_tau[], double *vpost, double *ipost);
void check_syn_type(SYN_TYPE types);

void update(int nstep, double *ic, neuron_t *cells, syn_t *syns, syn_t *bck_syns);

void add_isyn_bck(syn_t *syns);
void add_isyn(syn_t *syns);

void update_neurons(neuron_t *cells, int nstep);
void update_syns_no_delay(syn_t *syns, int *id_fired_pre);
void update_syns_delay(syn_t *syns, neuron_t *cells);
void gen_bck_spike(syn_t *bck_syns, int *id_fired_bck);

double *f_dv(double *v, void *arg_neuron, void *arg_null);
double *f_du(double *u, void *arg_neuron, void *arg_null);
double *f_dr_syns_no_delay(double *r, void *arg_syn, void *arg_fired);
double *f_dr_syns_delay(double *r, void *arg_syn, void *arg_cell);
double *f_dz_syns_delay(double *z, void *arg_syn, void *arg_null);
double *f_dz_syns_no_delay(double *z, void *arg_syn, void *arg_null);
double *solve_deq_using_euler(double* (*f) (double*, void*, void*), int N, double *x, void *arg1, void *arg2);
double *solve_deq_using_rk4(double* (*f) (double*, void*, void*), int N, double *x, void *arg1, void *arg2);
void append_spike(int nstep, int *num_spk, int **t_spk);
double get_avg(int num_x, double *x, int *is_target);
void reset_spike(neuron_t *cells);

// mathematical functions
void get_Kuramoto_order_params(int len, neuron_t *cells, int *is_target, double *rK, double *psiK);
void get_spike_phase(int n_spk, int nmax, int *nsteps, double *phase);
double *get_fft(int len, double *x_in);
double *get_fft_freq(int len, double fs);
void get_fft_summary(double *vm, double sample_rate, double t_range[2], arg_t *fft_res);
void downsampling(int len, double *y_org, double fs, arg_t *new_vars);
void get_summary(int max_step, double *vm, neuron_t *cells, int *targets, res_t *summary);

void free_summary(res_t *summary);
void free_neurons(neuron_t *cells);
void free_syns(syn_t *syns);
void free_rand_stream();
void destroy_mkl_buffers();

void init_network(network_info_t *info, neuron_t *cells, syn_t *syns, syn_t *bck_syns);
int *gen_types(int num, double *ratio);
double *linspace(double x0, double x1, int len_x);

#endif