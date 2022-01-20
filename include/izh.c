#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

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
#define malloc_c(size) MKL_malloc(size, 64)
#define calloc_c(len, var_size) MKL_calloc(len, var_size, 64)
#define free_c(arr) MKL_free(arr)

#else

#define malloc_c(size) malloc(size)
#define calloc_c(len, var_size) calloc(len, var_size)
#define free_c(arr) free(arr)

#endif

#define PI 3.14159265359
#define _block_size 500
VSLStreamStatePtr rand_stream;


// TODO: TWO cell interaction 확인해봐야될듯? inh가 잘 안들어가는 느낌
// TODO: 특정 initializing 함수들이 call됬는지 확인 변수 있으면 좋을듯
// TODO: compile bash도 따로 하나 general하게 만들면 좋을듯


double _dt = 0.005; // simulation time step
double _R  = 1;  // 
const double default_cell_params[4][4]= {
            {0.02, 0.2, -65, 8},    // RS
            {0.1, 0.2, -65, 2},   // FS
            {0.02, 0.2, -55, 4},    // IB
            {0.02, 0.2, -50, 2}};   // CH
const double default_syn_veq[4] = {0, -80, 0, 0};
const double default_syn_tau[4] = {5, 6, 5, 5};


void init_random_stream(long int seed)
{
    vslNewStream(&rand_stream, VSL_BRNG_MT19937, seed);
    init_genrand64(seed);
    // omp_set_num_threads(4);s
}

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
    cells->id_fired = (int*) malloc(sz_i * num_cells);

    for (int n=0; n<num_cells; n++){
        cells->t_fired[n] = (int*) malloc(sz_i * _block_size);
        cells->t_fired[n][0] = -1;
        cells->id_fired[n] = -1;
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
        syns->r = (double*) calloc_c(num_pres, sz_d);
        syns->inv_tau = (double*) malloc_c(sz_d * num_syns);
        syns->ptr_r = NULL;
    } else {
        syns->r = (double*) calloc_c(num_pres, sz_d);
        syns->inv_tau = (double*) malloc_c(sz_d * num_pres);
        syns->ptr_r = (double**) malloc(sizeof(double*) * num_syns);
    }

    syns->weight = (double*) malloc_c(sz_d * num_syns);
    syns->veq    = (double*) calloc_c(num_syns, sz_d);

    syns->ptr_vpost = (double**) malloc(sizeof(double*) * num_syns);
    syns->ptr_ipost = (double**) malloc(sizeof(double*) * num_syns);

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
            syns->weight[id_syn] = ntk->strength[n][i];

            int n_post = ntk->adj_list[n][i];
            syns->ptr_vpost[id_syn] = vpost + n_post;
            syns->ptr_ipost[id_syn] = ipost + n_post;

            syns->id_pre_neuron[id_syn]  = n;
            syns->id_post_neuron[id_syn] = n_post;

            if (syns->type == DELAY){
                syns->inv_tau[id_syn] = 1/cell_tau[pre_type];
            } else {
                syns->ptr_r[id_syn] = syns->r + n;
            }

            id_syn++;
        }

        if (syns->type != DELAY){
            syns->inv_tau[n] = 1/cell_tau[pre_type];
        }
    }
}


void update_no_delay(int nstep, double *ic, neuron_t *cells, syn_t *syns, syn_t *bck_syns)
{
    if (ic == NULL){
        memset(cells->ic, 0, sz_d*cells->num_cells);
    } else {
        memcpy(cells->ic, ic, sz_d*cells->num_cells);
    }

    int *id_fired_bck = (int*) malloc(sizeof(int) * (bck_syns->num_pres));
    gen_bck_spike(bck_syns, id_fired_bck);
    update_syns_no_delay(bck_syns, id_fired_bck);
    add_isyn_bck(bck_syns);
    free(id_fired_bck);

    update_syns_no_delay(syns, cells->id_fired);
    add_isyn(syns);

    update_neurons(cells, nstep);
}


// void update_delay()
// {

// }


void add_isyn_bck(syn_t *syns)
{
    /*** ic -= weight * r * vpost ***/
    int N = syns->num_syns;

    double *rsyns = (double*) malloc_c(sz_d * N);
    read_ptr(N, rsyns, syns->ptr_r);

    double *isyn = (double*) malloc_c(sz_d * N);
    read_ptr(N, isyn, syns->ptr_vpost);

    vdMul(N, rsyns, isyn, isyn); // W*r*v
    vdMul(N, syns->weight, isyn, isyn);

    for (int n=0; n<N; n++){
        *(syns->ptr_ipost[n]) -= isyn[n];
    }

    free_c(isyn);
    free_c(rsyns);
}


void add_isyn(syn_t *syns)
{
    /*** ic -= weight * r * (vpost - veq) ***/
    int N = syns->num_syns;

    // int nid = 0;
    // while (syns->id_pre_neuron[nid] < 80){
    //     nid ++;
    // }

    double *rsyns = (double*) malloc_c(sz_d * N);
    read_ptr(N, rsyns, syns->ptr_r);

    double *isyn = (double*) malloc_c(sz_d * N);
    read_ptr(N, isyn, syns->ptr_vpost);
    cblas_daxpy(N, -1, syns->veq, 1, isyn, 1);  // vpost - veq
    vdMul(N, syns->weight, isyn, isyn);
    vdMul(N, rsyns, isyn, isyn);

    for (int n=0; n<N; n++){
        *(syns->ptr_ipost[n]) -= isyn[n];
    }

    free_c(isyn);
    free_c(rsyns);
}


void add_isyn_delay(syn_t *syns)
{
    /*** ic -= weight * (vpost - veq) ***/
    int N = syns->num_syns;

    double *isyn = (double*) malloc_c(sz_d * N);
    read_ptr(N, isyn, syns->ptr_vpost);

    cblas_daxpy(N, -1, syns->veq, 1, isyn, 1);  // vpost - veq
    vdMul(N, syns->weight, isyn, isyn);
    vdMul(N, syns->r, isyn, isyn);

    for (int n=0; n<N; n++){
        *(syns->ptr_ipost)[n] -= isyn[n];
    }

    free_c(isyn);
}


void update_neurons(neuron_t *cells, int nstep)
{
    int N = cells->num_cells;
    double *dv = solve_deq_using_rk4(f_dv, N, cells->v, (void*) cells, NULL);
    double *du = solve_deq_using_rk4(f_du, N, cells->u, (void*) cells, NULL);

    cblas_daxpy(N, 1, dv, 1, cells->v, 1);
    cblas_daxpy(N, 1, du, 1, cells->u, 1);

    // check spike
    int *ptr_id = cells->id_fired;
    for (int n=0; n<N; n++){
        if (cells->v[n] > 30){
            cells->v[n] = cells->c[n];
            cells->u[n] += cells->d[n];

            // append spike
            *ptr_id++ = n;
            append_spike(nstep, cells->num_spk+n, cells->t_fired+n);
        }
    }
    *ptr_id = -1;

    free_c(dv);
    free_c(du);
}


void update_syns_no_delay(syn_t *syns, int *id_fired_pre)
{   
    double *dr = solve_deq_using_rk4(f_dr_syns_no_delay,
                        syns->num_pres, syns->r, (void*) syns, (void*) id_fired_pre);
    cblas_daxpy(syns->num_pres, 1, dr, 1, syns->r, 1);
    free_c(dr);
}


void update_syns_delay(syn_t *syns)
{
    /*** TODO: fill this function ***/
}


void gen_bck_spike(syn_t *bck_syns, int *id_fired_bck)
{
    int N = bck_syns->num_pres;
    double *p = (double*) malloc_c(sz_d * N);
    
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rand_stream, N, p, 0, 1);
    
    int *ptr_id = id_fired_bck;
    for (int n=0; n<N; n++){
        if (p[n] < bck_syns->p_fire[n]){
            *ptr_id++ = n;
        }
    }
    *ptr_id = -1;

    free_c(p);
}


double *f_dv(double *v, void *arg_neuron, void *arg_null)
{
    /*** dv/dt = 0.04*v^2 + 5*v - u + I + 140 ***/
    neuron_t *cells = (neuron_t*) arg_neuron;
    int N = cells->num_cells;

    double *dv = (double*) malloc_c(sz_d * N);
    
    vdMul(N, v, v, dv);
    cblas_daxpby(N, 5*_dt, v, 1, 0.04*_dt, dv, 1);
    cblas_daxpy(N, -_dt, cells->u, 1, dv, 1);
    cblas_daxpy(N, _dt, cells->ic, 1, dv, 1);

    double c=140 *_dt;
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

    int *ptr_id = (int*) arg_fired;
    while (*ptr_id > -1){
        dr[*ptr_id++] += _R; // dx = deltafn(=_R/_dt) * _dt
    }

    vdMul(syns->num_pres, syns->inv_tau, dr, dr);

    return dr;
}


// double *f_dr_syns_delay(double *r, void *arg_syn, void *arg_fired)
// {
//     /*** TODO: fill the function***/
//     syn_t *syns = (syn_t*) arg_syn;
//     int *id_fired = (int *) arg_fired;
// }


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

    cblas_daxpby(N, 1./3, dx2, 1, 1./6, dx, 1);
    cblas_daxpy(N,  1./3, dx3, 1, dx, 1);
    cblas_daxpy(N,  1./3, dx4, 1, dx, 1);

    free_c(dx2);
    free_c(dx3);
    free_c(dx4);

    return dx;
}


void get_Kuramoto_order_params(int len, neuron_t *cells, int *is_target, double *rK, double *psiK)
{   
    int num_cells=cells->num_cells;

    double *sum_real = (double*) calloc(len, sz_d);
    double *sum_imag = (double*) calloc(len, sz_d);
    int N = 0;

    for (int n=0; n<num_cells; n++){
        if ((is_target != NULL) && (is_target[n] == 0)){
            continue;
        }
        N++;

        double *phase = (double*) malloc(sz_d * len);
        get_spike_phase(cells->num_spk[n], len, cells->t_fired[n], phase);

        double *tmp_real = (double*) malloc(sz_d * len);
        double *tmp_imag = (double*) malloc(sz_d * len);

        vdCos(len, phase, tmp_real);
        vdSin(len, phase, tmp_imag);
        // for (int n=0; n<len; n++){
        //     tmp_real[n] = cos(phase[n]);
        //     tmp_imag[n] = sin(phase[n]);

        //     // sum_real[n] += tmp_real
        // }

        vdAdd(len, tmp_real, sum_real, sum_real);
        vdAdd(len, tmp_imag, sum_imag, sum_imag);
        
        free(phase);
        free(tmp_real);
        free(tmp_imag);
    }

    cblas_dscal(len, 1./N, sum_real, 1);
    cblas_dscal(len, 1./N, sum_imag, 1);

    // get Kuramoto order parameter
    for (int i=0; i<len; i++){
        psiK[i] = atan2(sum_real[i], sum_imag[i]); // [-phi, phi]
        rK[i] = sqrt(sum_real[i]*sum_real[i] + sum_imag[i]*sum_imag[i]);
    }

    free(sum_real);
    free(sum_imag);
}


void get_spike_phase(int n_spk, int nmax, int *nsteps, double *phase)
{
    int i, *n0, *n1, tmp=0;

    memset(phase, 0, sizeof(double) * nmax);

    if (n_spk == 0){
        return;
    }

    n0 = &tmp;
    n1 = nsteps;
    n_spk--;

    for (i=0; i<nmax; i++){
        if (i == *n1) {
            if (n_spk == 0){
                return;
            }
            n0 = n1;
            n1++;
        }
        phase[i] = 2*PI * (i - *n0)/(*n1 - *n0);
    }
}

double get_avg(int num_x, double *x, int *is_target)
{
    double x_avg = 0;

    int N = 0;
    for (int n=0; n<num_x; n++){
        if ((is_target != NULL) && (is_target[n] == 0)){
            continue;
        }
        x_avg += x[n];
        N++;
    }
    x_avg /= N;
    return x_avg;
}


void reset_spike(neuron_t *cells)
{
    for (int n=0; n<cells->num_cells; n++){
        for (int i=0; i<cells->num_spk[n]; i++){
            cells->t_fired[n][i] = 0;
        }
        cells->num_spk[n] = 0;
    }
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


void read_ptr(int num_x, double *x, double **ptr_x)
{
    // #pragma omp parallel for if (num_x > 1000)
    for (int n=0; n<num_x; n++){
        x[n] = *(ptr_x[n]);
    }
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
    free(cells->id_fired);
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

    free(syns->id_pre_neuron);
    free(syns->id_post_neuron);

    free(syns->ptr_r);
    free(syns->ptr_vpost);
    free(syns->ptr_ipost);
    
    if (syns->type == BACKGROUND){
        free(syns->p_fire);
    } else if (syns->type == DELAY) {
        free(syns->delay);
    }
}


void free_rand_stream()
{
    vslDeleteStream(&rand_stream);
}