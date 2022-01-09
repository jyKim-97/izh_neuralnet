#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mkl.h>
#include <mkl_vsl.h>
#include "mt64.h"
#include "izh.h"

// set simulation env variable

#define USE_MKL_MEM
#define sz_d sizeof(double)
#define sz_i sizeof(int)

#ifdef USE_MKL_MEM

#define ALIGN_SIZE
#define malloc_c(size) MKL_malloc(size, 32)
#define calloc_c(len, var_size) MKL_calloc(len, var_size, 32)
#define free_c(arr) MKL_free(arr)

#else

#define malloc_c(size) malloc(size)
#define calloc_c(len, var_size) calloc(len, var_size)
#define free_c(arr) free(arr)

#endif

#define PI 3.14159265359
#define _block_size 500
VSLStreamStatePtr rand_stream;

double _dt = 0.005; // simulation time step
double _R  = 0.001;  // 

void init_random_stream(long int seed)
{
    vslNewStream(&rand_stream, VSL_BRNG_MT19937, seed);
    init_genrand64(seed);
}

// Background synapse firing probability 처리?

void init_cell_vars(neuron_t *cells, int num_cells, double cell_params[][4], int *cell_types)
{
    cells->num_cells = num_cells;
    
    cells->v  = (double*) malloc_c(sz_d * num_cells);
    cells->u  = (double*) malloc_c(sz_d * num_cells);
    cells->ic = (double*) malloc_c(sz_d * num_cells);
    cells->a  = (double*) malloc_c(sz_d * num_cells);
    cells->b  = (double*) malloc_c(sz_d * num_cells);
    cells->c  = (double*) malloc_c(sz_d * num_cells);
    cells->d  = (double*) malloc_c(sz_d * num_cells);

    // allocate cell params
    for (int n=0; n<num_cells; n++){
        cells->v[n] = -70;
        cells->u[n] = 0;
        // set parameters
        int ctp = cell_types[n];
        cells->a[n] = cell_params[ctp][0];
        cells->b[n] = cell_params[ctp][1];
        cells->c[n] = cell_params[ctp][2];
        cells->d[n] = cell_params[ctp][3];
    }

    // init spike parameters
    cells->num_spk = (int*) calloc(num_cells, sz_i);
    cells->t_fired = (int**) malloc(sizeof(int*) * num_cells);
    cells->id_fired_neuron = (int*) malloc(sz_i * num_cells);

    for (int n=0; n<num_cells; n++){
        cells->t_fired[n] = (int*) malloc(sz_i * _block_size);
        cells->t_fired[n][0] = -1;
        cells->id_fired_neuron[n] = -1;
    }
}


void init_syn_vars(syn_t *syns, int num_pres, SYN_TYPE type, ntk_t *ntk, double cell_veq[], double cell_tau[], double *vpost, double *ipost)
{
    int num_syns = 0;
    for (int n=0; n<num_pres; n++){
        num_syns += ntk->num_edges[n];
    }
    syns->num_syns = num_syns;
    syns->num_pres = num_pres;
    syns->type = type;

    // create objs
    if (syns->type == DELAY){
        syns->r = (double*) mallo_c(sz_d * num_syns);
        syns->inv_tau = (double*) mallo_c(sz_d * num_syns);
        syns->ptr_r = NULL;
    } else {
        syns->r = (double*) malloc_c(sz_d * num_pres);
        syns->inv_tau = (double*) mallo_c(sz_d * num_pres);
        syns->ptr_r = (double**) malloc(sizeof(double*) * num_syns);
    }

    syns->weight = (double*) malloc_c(sz_d * num_syns);
    syns->veq    = (double*) malloc_c(sz_d * num_syns);

    syns->ptr_vpost = (double*) malloc(sizeof(double*) * num_syns);
    syns->ptr_ipost = (double*) malloc(sizeof(double*) * num_syns);

    syns->id_pre_neuron  = (int*) malloc(sz_i * num_syns);
    syns->id_post_neuron = (int*) malloc(sz_i * num_syns);

    if (syns->type == BACKGROUND){
        syns->p_fire = (double*) malloc(sz_d * num_pres);
        syns->delay = NULL;
    } else if (syns->type == DELAY){
        syns->delay = (int*) malloc(sz_i * num_syns);
        syns->p_fire = NULL;
    } else {
        syns->p_fire = NULL;
        syns->delay  = NULL;
    }

    // set varaibles
    int id_syn = 0;
    for (int n=0; n<num_pres; n++){
        int pre_type = ntk->node_types[n];
        for (int i=0; i<ntk->num_edges[n]; i++){
            syns->veq[id_syn] = cell_veq[pre_type];
            // syns->inv_tau[id_syn] = 1/cell_tau[pre_type];
            syns->weight[id_syn] = ntk->strength[n][i];

            int n_post = ntk->adj_list[n][i];
            syns->ptr_vpost[id_syn] = vpost + n_post;
            syns->ptr_ipost[id_syn] = vpost + n_post;

            syns->id_pre_neuron[id_syn]  = n;
            syns->id_post_neuron[id_syn] = n;

            if (syns->type == DELAY){
                syns->inv_tau[id_syn] = 1/cell_tau[pre_type];
            } else {
                syns->ptr_r[id_syn] = syns->r + n;
            }

            id_syn++;
        }

        if (syns->type < 2){
            syns->inv_tau[n] = 1/cell_veq[pre_type];
        }
    }
}


void add_isyn_bck(syn_t *syns)
{
    /*** ic -= weight * (vpost - veq) ***/
    int N = syns->num_syns;

    // double *vpost = read_ptr(N, syns->ptr_vpost);
    double *rsyns = read_ptr(N, syns->ptr_r);

    double *isyn = read_ptr(N, syns->ptr_vpost); // vpost
    // cblas_daxpy(N, -1, syns->veq, 1, isyn, 1);  // vpost - veq
    vdMul(N, syns->weight, isyn, isyn);

    for (int n=0; n<N; n++){
        (*syns->ptr_ipost)[n] -= isyn[n];
    }

    free_c(isyn);
    free_c(rsyns);
}


void add_isyn(syn_t *syns)
{
    /*** ic -= weight * (vpost - veq) ***/
    int N = syns->num_syns;

    // double *vpost = read_ptr(N, syns->ptr_vpost);
    double *rsyns = read_ptr(N, syns->ptr_r);

    double *isyn = read_ptr(N, syns->ptr_vpost); // vpost
    cblas_daxpy(N, -1, syns->veq, 1, isyn, 1);  // vpost - veq
    vdMul(N, syns->weight, isyn, isyn);
    vdMul(N, rsyns, isyn, isyn);

    for (int n=0; n<N; n++){
        (*syns->ptr_ipost)[n] -= isyn[n];
    }

    free_c(isyn);
    free_c(rsyns);
}


void add_isyn_delay(syn_t *syns)
{
    /*** ic -= weight * (vpost - veq) ***/
    int N = syns->num_syns;

    double *isyn = read_ptr(N, syns->ptr_vpost); // vpost
    cblas_daxpy(N, -1, syns->veq, 1, isyn, 1);  // vpost - veq
    vdMul(N, syns->weight, isyn, isyn);
    vdMul(N, syns->r, isyn, isyn);

    for (int n=0; n<N; n++){
        (*syns->ptr_ipost)[n] -= isyn[n];
    }

    free_c(isyn);
}


void update_neurons(neuron_t *cells, int nstep, int *id_fired_neuron)
{
    int N = cells->num_cells;
    double *dv = solve_deq_using_rk4(f_dv, N, cells->v, (void*) cells, NULL);
    double *du = solve_deq_using_rk4(f_dv, N, cells->v, (void*) cells, NULL);

    cblas_daxpy(N, 1, dv, 1, cells->v, 1);
    cblas_daxpy(N, 1, dv, 1, cells->v, 1);

    // check spike
    double *ptr_v = cells->v;
    int *ptr_id = id_fired_neuron;
    for (int n=0; n<N; n++){
        if (*ptr_v > 30){
            *ptr_id++ = n;
            *ptr_v = cells->c[n];
            cells->u[n] += cells->d[n];

            // append spike
            append_spike(nstep, cells->num_spk+n, cells->t_fired+n);
        }
    }
    *ptr_id = -1;

    free_c(dv);
    free_c(du);
}


void update_syns_no_delay(syn_t *syns, int *id_fired_neuron)
{
    double *dr = solve_deq_using_rk4(f_dr_syns_no_delay,
                        syns->num_pres, syns->r, (void*) syns, (void*) id_fired_neuron);
    cblas_daxpy(syns->num_pres, 1, dr, 1, syns->r, 1);
    free_c(dr);
}


void update_syns_delay(syn_t *syns)
{
    /*** TODO: fill this function ***/
}


int *gen_bck_spike(syn_t *bck_syns)
{
    int N = bck_syns->num_pres;
    double *p = (double*) malloc_c(sz_d * N);
    
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rand_stream, N, p, 0, 1);
    
    int *id_fired_bck = (int*) calloc_c(N, sz_i);
    for (int n=0; n<N; n++){
        if (p[n] < bck_syns->p_fire[n]){
            id_fired_bck[n] = 1;
        }
    }

    free_c(p);
    return id_fired_bck;
}


double *f_dv(double *v, void *arg_neuron, void *arg_null)
{
    /*** dv/dt = 0.04*v^2 + 5*v - u + I + 140 ***/
    neuron_t *cells = (neuron_t*) arg_neuron;
    int N = cells->num_cells;

    double *dv = (double*) malloc_c(sz_d * N);
    
    vdPow(N, v, 2, dv);
    cblas_daxpby(N, 5*_dt, v, 1, 0.04*_dt, dv, 1);
    cblas_daxpy(N, -_dt, cells->u, 1, dv, 1);
    cblas_daxpy(N, _dt, cells->ic, 1, dv, 1);

    double c=140*_dt;
    for (int n=0; n<N; n++){
        dv[n] += c;
    }

    return dv;
}


double *f_du(double *u, void *arg_neuron, void *arg_null)
{
    /*** du/dt = a*(b*v-u) ***/
    neuron_t *cells = (neuron_t*) arg_neuron;
    int N = cells->num_cells;

    double *du = (double*) malloc_c(sz_d * N);

    vdMul(N, cells->b, cells->v, du);
    cblas_daxpy(N, -1, u, 1, du, 1);
    vdMul(N, cells->a, du, du);
    cblas_dscal(N, _dt, du, 1);

    return du;
}


double *f_dr_syns_no_delay(double *r, void *arg_syn, void *arg_fired)
{
    /*** dr/dt = (-r + R\delta) / tau ***/
    syn_t *syns = (syn_t*) arg_syn;
    double *dr = (double*) malloc_c(sz_d * syns->num_pres);

    memcpy(dr, r, sz_d*syns->num_pres);
    cblas_dscal(syns->num_pres, -_dt, dr, 1);

    int *ptr_id = arg_fired;
    while (*ptr_id > 0){
        dr[*ptr_id++] += _R; // _R/_dt *
    }

    vdMul(syns->num_pres, syns->inv_tau, dr, dr);

    return dr;
}


double *f_dr_syns_delay(double *r, void *arg_syn, void *arg_fired)
{
    /*** TODO: fill the function***/
    syn_t *syns = (syn_t*) arg_syn;
    int *id_fired_neuron = (int *) arg_fired;
}


double *solve_deq_using_rk4(double* (*f) (double*, void*, void*), int N, double *x, void *arg1, void *arg2)
{
    // solve differential equation usign Runge-Kutta 4th order method

    double *dx = f(x, arg1, arg2);

    double *xtmp = (double*) malloc_c(sz_d * N);
    vdAdd(N, x, dx, xtmp);
    double *dx2 = f(xtmp, arg1, arg2);

    vdAdd(N, x, dx2, xtmp);
    double *dx3 = f(xtmp, arg1, arg2);

    vdAdd(N, x, dx3, xtmp);
    double *dx4 = f(xtmp, arg1, arg2);
    free_c(xtmp);

    cblas_daxpby(N, 1/3, dx2, 1, 1/6, dx, 1);
    cblas_daxpy(N,  1/3, dx3, 1, dx, 1);
    cblas_daxpy(N,  1/3, dx4, 1, dx, 1);

    free_c(dx2);
    free_c(dx3);
    free_c(dx4);

    return dx;
}


void get_avg(int num_x, int **ptr_x)
{
    double x = 0;
    for (int n=0; n<num_x; n++){
        x += **ptr_x++;
    }
    return x/num_x;
}


void append_spike(int nstep, int *num_spk, int **t_spk)
{
    (*t_spk)[*num_spk] = nstep;
    (*num_spk)++;
    if (((*num_spk) % _block_size == 0) && (*num_spk > 0)){
        *t_spk = (int*) realloc(*t_spk, (*num_spk + _block_size) * sz_d);
    }
    (*t_spk)[*num_spk] = -1;
}


double *read_ptr(int num_x, double **ptr_x)
{
    double *x = (double*) malloc_c(sz_d * num_x);

    for (int n=0; n<num_x; n++){
        x[n] = **ptr_x++;
    }

    return x;
}


void free_neurons(neuron_t *cells)
{
    free_c(cells->v);
    free_c(cells->u);
    free_c(cells->ic);
    free_c(cells->a);
    free_c(cells->b);
    free_c(cells->c);
    free_c(cells->d);
    free(cells->num_spk);
    free(cells->id_fired_neuron);
    for (int n=0; n<cells->num_cells; n++){
        free(cells->t_fired[n]);
    }
    free(cells->t_fired);
}


void free_syns(syn_t *syns)
{
    free_c(syns->r);
    free_c(syns->veq);
    free_c(syns->weight);
    free_c(syns->inv_tau);

    free(syns->ptr_r);
    free(syns->ptr_vpost);
    free(syns->ptr_ipost);
    
    if (syns->type == BACKGROUND){
        free(syns->p_fire);
    } else if (syns->type == DELAY) {
        free(syns->delay);
    }
}