#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <stdarg.h>
#include <stdbool.h>
#include <omp.h>

#include "mkl.h"
#include "mkl_vsl.h"
#include "mkl_dfti.h"
#include "mt64.h"
#include "izh.h"

// set simulation env variable
// #define USE_MKL_MEM
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
#define _block_size 2000
VSLStreamStatePtr rand_stream;


double _dt = 0.005; // simulation time step
const double default_cell_params[MAX_TYPE][4]= {
        {0.02, 0.2, -65, 8},    // RS
            {0.1, 0.2, -65, 2},   // FS
            {0.02, 0.2, -55, 4},    // IB
            {0.02, 0.2, -50, 2},
            {0, 0, 0, 0},
            {0, 0, 0, 0}};   // CH
const double default_syn_veq[MAX_TYPE] = {0, -80, 0, 0};
const double default_syn_tau[MAX_TYPE] = {5, 6, 5, 5};
double t_skip_phs = 5; // Using when calculate spike phase, ms
double fs = 2000;
double f_cut[2] = {4, 300};
bool pass_spk_id = true;


void init_random_stream(long int seed)
{
    vslNewStream(&rand_stream, VSL_BRNG_MT19937, seed);
    init_genrand64(seed);
}

void init_cell_vars(neuron_t *cells, int num_cells, double cell_params[MAX_TYPE][4], int *cell_types)
{
    cells->num_cells = num_cells;

    cells->v  = (double*) malloc_c(sz_d * num_cells);
    cells->u  = (double*) malloc_c(sz_d * num_cells);
    // cells->ic = (double*) malloc_c(sz_d * num_cells);
    cells->ic = (double*) calloc_c(num_cells, sz_d);
    cells->a  = (double*) malloc_c(sz_d * num_cells);
    cells->b  = (double*) malloc_c(sz_d * num_cells);
    cells->c  = (double*) malloc_c(sz_d * num_cells);
    cells->d  = (double*) malloc_c(sz_d * num_cells);

    // allocate cell params
    for (int n=0; n<num_cells; n++){
        cells->v[n] = genrand64_normal(-70, 5);
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
    cells->id_fired = (int*) malloc(sz_i * (num_cells+1));

    for (int n=0; n<num_cells; n++){
        cells->t_fired[n] = (int*) malloc(sz_i * _block_size);
        cells->t_fired[n][0] = -1;
    }
    cells->id_fired[0] = -1;

    cells->types = (int*) malloc(sizeof(int) * num_cells);
    memcpy(cells->types, cell_types, sizeof(int) * num_cells);
}


void init_syn_vars(syn_t *syns, int num_pres, SYN_TYPE types, ntk_t *ntk, double cell_veq[], double cell_tau[], double *vpost, double *ipost)
{
    int num_syns = 0;
    for (int n=0; n<num_pres; n++){
        num_syns += ntk->num_edges[n];
    }
    syns->num_syns = num_syns;
    syns->num_pres = num_pres;

    check_syn_type(types);
    syns->type = types;

    // create objs
    if (syns->type & DELAY){
        syns->r = (double*) calloc_c(num_syns, sz_d);
        syns->inv_tau = (double*) malloc_c(sz_d * num_syns);
        syns->delay = (int*) malloc(sz_i * num_syns);
        syns->id_exp = (int*) calloc(num_syns, sz_i);
        syns->p_fire = NULL;

        for (int n=0; n<num_syns; n++) syns->delay[n] = -1;

        if (syns->type & STD){
            syns->tau_r = 100;
            syns->x = (double*) malloc_c(sz_d * num_syns);
            syns->z = (double*) calloc_c(num_syns, sz_d);
            for (int i=0; i<num_syns; i++) syns->x[i] = 1;
        }
    } else {
        syns->r = (double*) calloc_c(num_pres, sz_d);
        syns->inv_tau = (double*) malloc_c(sz_d * num_pres);
        syns->delay  = NULL;
        syns->id_exp = NULL;
        if (syns->type & NO_DELAY){
            syns->p_fire = NULL;
            if (syns->type & STD){
                syns->tau_r = 100;
                syns->x = (double*) malloc_c(sz_d * num_pres);
                syns->z = (double*) calloc_c(num_pres, sz_d);
                for (int i=0; i<num_pres; i++) { syns->x[i] = 1; }
            }
        } else { // background
            syns->p_fire = (double*) malloc(sz_d * num_pres);
            syns->delay = NULL;
        }
    }

    syns->veq = (double*) calloc_c(num_pres, sz_d);
    syns->weight = (double*) malloc_c(sz_d * num_syns);
    syns->id_pre_neuron = (int*) malloc(sz_i * num_syns);
    syns->id_post_neuron = (int*) malloc(sz_i * num_syns);
    syns->ptr_ipost = (double**) malloc(sizeof(double*) * num_syns);
    syns->ptr_vpost = (double**) malloc(sizeof(double*) * num_syns);

    // set varaibles    
    int id_syn = 0;
    for (int n_pre=0; n_pre<num_pres; n_pre++){

        int pre_type = ntk->node_types[n_pre];
        syns->veq[n_pre] = cell_veq[pre_type];
        if (syns->type ^ DELAY) syns->inv_tau[n_pre] = 1/cell_tau[pre_type];

        for (int i=0; i<ntk->num_edges[n_pre]; i++){
            syns->weight[id_syn] = ntk->strength[n_pre][i];

            int n_post = ntk->adj_list[n_pre][i];
            syns->ptr_ipost[id_syn] = ipost + n_post;
            syns->ptr_vpost[id_syn] = vpost + n_post;
            syns->id_pre_neuron[id_syn]  = n_pre;
            syns->id_post_neuron[id_syn] = n_post;

            if (syns->type & DELAY) syns->inv_tau[id_syn] = 1/cell_tau[pre_type];
            id_syn++;
        }
    }
}


void check_syn_type(SYN_TYPE types){
    for (int n=0; n<3; n++){
        if (types & 0x01<<n){
            return;
        }
    }
    printf("Synapse type is not specified! check type = %d\n", types);
}


void update(int nstep, double *ic, neuron_t *cells, syn_t *syns, syn_t *bck_syns)
{
    if (ic==NULL){
        memset(cells->ic, 0, sz_d*cells->num_cells);
    } else {
        memcpy(cells->ic, ic, sz_d*cells->num_cells);
    }

    int *id_fired_bck = (int*) malloc(sizeof(int) * (bck_syns->num_pres+1));
    gen_bck_spike(bck_syns, id_fired_bck);
    update_syns_no_delay(bck_syns, id_fired_bck);
    add_isyn_bck(bck_syns);
    free(id_fired_bck);

    if (syns->type & DELAY){
        update_syns_delay(syns, cells);
    } else {
        update_syns_no_delay(syns, cells->id_fired);
    }
    add_isyn(syns);

    update_neurons(cells, nstep);
}


void add_isyn_bck(syn_t *syns)
{
    /*** ic -= weight * r * vpost ***/
    for (int n=0; n<syns->num_syns; n++){
        int id = syns->id_pre_neuron[n];
        *(syns->ptr_ipost[n]) -= syns->weight[n] * syns->r[id] * *(syns->ptr_vpost[n]);
    }
}


void add_isyn(syn_t *syns)
{
    if (syns->type & NO_DELAY){
        for (int n=0; n<syns->num_syns; n++){
            int id = syns->id_pre_neuron[n];
            double rdv = syns->r[id] * (*(syns->ptr_vpost[n]) - syns->veq[id]);
            *(syns->ptr_ipost[n]) -= syns->weight[n]*rdv;
        }

    } else if (syns->type & DELAY){
        for (int n=0; n<syns->num_syns; n++){
            int id = syns->id_pre_neuron[n];
            double dv = *(syns->ptr_vpost[n]) - syns->veq[id];
            *(syns->ptr_ipost[n]) -= syns->weight[n]*syns->r[n]*dv;
        }
    }
}


void update_neurons(neuron_t *cells, int nstep)
{
    int N = cells->num_cells;
    cells->nstep = nstep;
    
    double *dv = solve_deq_using_rk4(f_dv, N, cells->v, (void*) cells, NULL);
    double *du = solve_deq_using_rk4(f_du, N, cells->u, (void*) cells, NULL);
    // double *dv = solve_deq_using_rk2(f_dv, N, cells->v, (void*) cells, NULL);
    // double *du = solve_deq_using_rk2(f_du, N, cells->u, (void*) cells, NULL);
    // double *dv = solve_deq_using_euler(f_dv, N, cells->v, (void*) cells, NULL);
    // double *du = solve_deq_using_euler(f_du, N, cells->u, (void*) cells, NULL);

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
    int N = syns->num_pres;
    // double *dr = solve_deq_using_euler(f_dr_syns_no_delay, N, syns->r, (void*) syns, (void*) id_fired_pre);
    // double *dr = solve_deq_using_rk2(f_dr_syns_no_delay, N, syns->r, (void*) syns, (void*) id_fired_pre);
    double *dr = solve_deq_using_rk4(f_dr_syns_no_delay, N, syns->r, (void*) syns, (void*) id_fired_pre);
    cblas_daxpy(N, 1, dr, 1, syns->r, 1);

    if (syns->type & STD){
        double *dz = solve_deq_using_euler(f_dz_syns_no_delay, N, syns->z, (void*) syns, NULL);
        cblas_daxpy(N, 1, dz, 1, syns->z, 1);
        cblas_daxpy(N, 1, dr, 1, dz, 1); // dz = dr + dz, dx = -dr
        cblas_daxpy(N, -1, dz, 1, syns->x, 1); // dx + dr + dz = 0
        free_c(dz);
    }

    free_c(dr);
}


void update_syns_delay(syn_t *syns, neuron_t *cells)
{
    int N = syns->num_syns;
    // double *dr = solve_deq_using_euler(f_dr_syns_delay, N, syns->r, (void*) syns, (void*) cells);
    // double *dr = solve_deq_using_rk2(f_dr_syns_delay, N, syns->r, (void*) syns, (void*) cells);
    double *dr = solve_deq_using_rk4(f_dr_syns_delay, N, syns->r, (void*) syns, (void*) cells);
    cblas_daxpy(N, 1, dr, 1, syns->r, 1);

    if (syns->type & STD){
        double *dz = solve_deq_using_euler(f_dz_syns_delay, N, syns->z, (void*) syns, NULL);
        cblas_daxpy(N, 1, dz, 1, syns->z, 1);
        cblas_daxpy(N, 1, dr, 1, dz, 1); // dr = dr + dz, dx = -dr
        cblas_daxpy(N, -1, dz, 1, syns->x, 1); // dx + dr + dz = 0
        free_c(dz);
    }

    free_c(dr);
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


double *f_dz_syns_delay(double *z, void *arg_syn, void *arg_null)
{
    syn_t *syns = (syn_t*) arg_syn;
    int N = syns->num_syns;
    double *dz = (double*) malloc_c(sz_d*N);

    memcpy(dz, z, sz_d*N);

    double *dy = (double*) malloc_c(sz_d * N);
    vdMul(N, syns->inv_tau, syns->r, dy);
    cblas_daxpby(N, _dt, dy, 1, -_dt/syns->tau_r, dz, 1);

    free(dy);
    return dz;
}


double *f_dz_syns_no_delay(double *z, void *arg_syn, void *arg_null)
{
    syn_t *syns = (syn_t*) arg_syn;

    int N = syns->num_pres;
    double *dz = (double*) malloc_c(sz_d * N);

    memcpy(dz, z, sz_d*N);

    double *dy = (double*) malloc_c(sz_d * N);
    vdMul(N, syns->inv_tau, syns->r, dy);

    cblas_daxpby(N, _dt, dy, 1, -_dt/syns->tau_r, dz, 1);
 
    free(dy);

    return dz;
}


double *f_dr_syns_no_delay(double *r, void *arg_syn, void *arg_fired)
{
    /*** dr/dt = (-r + R\delta) / tau ***/
    syn_t *syns = (syn_t*) arg_syn;
    double *dr = (double*) malloc(sz_d * syns->num_pres);

    memcpy(dr, r, sz_d*syns->num_pres);
    cblas_dscal(syns->num_pres, -_dt, dr, 1);

    int *ptr_id = (int*) arg_fired;
    while (*ptr_id > -1){
        if (syns->type & STD){
            dr[*ptr_id] += syns->x[*ptr_id]; // dx = deltafn(=1/_dt) * _dt
        } else {
            dr[*ptr_id] += 1; // dx = deltafn(=1/_dt) * _dt
        }

        ptr_id++;
    }

    double *dr_new = (double*) malloc_c(sz_d * syns->num_pres);
    vdMul(syns->num_pres, syns->inv_tau, dr, dr_new);
    free(dr);

    return dr_new;
}


double *f_dr_syns_delay(double *r, void *arg_syn, void *arg_cell)
{
    syn_t *syns = (syn_t*) arg_syn;
    neuron_t *cells = (neuron_t*) arg_cell;
    double *dr = (double*) malloc_c(sz_d * syns->num_syns);

    memcpy(dr, r, sz_d*syns->num_syns);
    cblas_dscal(syns->num_syns, -_dt, dr, 1);

    int nstep=cells->nstep;
    for (int n=0; n<syns->num_syns; n++){
        int id_pre = syns->id_pre_neuron[n];
        int delay = syns->delay[n];
        int n_spk = cells->num_spk[id_pre];
        int n_exp = syns->id_exp[n];

        if (n_exp < n_spk){
            int dn = nstep - delay - cells->t_fired[id_pre][n_exp];
            if (dn == 0){
                if (syns->type & STD){
                    dr[n] += syns->x[n];
                    // dr[n] += 1;
                } else {
                    dr[n] += 1;
                }

                if (pass_spk_id){
                    syns->id_exp[n]++;
                }
            }
        }
    }

    // if (pass_spk_id)

    double *dr_new = (double*) malloc_c(sz_d * syns->num_syns);
    vdMul(syns->num_syns, syns->inv_tau, dr, dr_new);
    free(dr);

    return dr_new;
}


double *solve_deq_using_euler(double* (*f) (double*, void*, void*), int N, double *x, void *arg1, void *arg2)
{
    double *dx = f(x, arg1, arg2);
    return dx;
}


double *solve_deq_using_rk4(double* (*f) (double*, void*, void*), int N, double *x, void *arg1, void *arg2)
{
    // solve differential equation usign Runge-Kutta 4th order method
    pass_spk_id = true;
    double *xtmp = (double*) malloc_c(sz_d * N);

    // calculate dx1
    double *dx = f(x, arg1, arg2);

    // dx2
    cblas_dscal(N, 1./2., dx, 1); // dx1' <- dx1/2
    vdAdd(N, x, dx, xtmp); // xtmp <- x + dx1/2 = x + dx1'
    pass_spk_id = false;
    double *dx2 = f(xtmp, arg1, arg2);

    // dx3
    cblas_dscal(N, 1./2., dx2, 1); // dx2' <- dx2/2
    vdAdd(N, x, dx2, xtmp); // xtmp <- x + dx2/2 = x + dx2'
    double *dx3 = f(xtmp, arg1, arg2);

    // dx4
    vdAdd(N, x, dx3, xtmp); // xtmp <- x + dx3
    double *dx4 = f(xtmp, arg1, arg2);
    free_c(xtmp);

    // dx = (dx1 + 2dx2 + 2dx3 + dx4)/6 = (2dx1' + 4dx2' + 2dx3 + dx4)/6
    cblas_daxpby(N, 2./3., dx2, 1, 1./3., dx, 1);
    cblas_daxpy(N,  1./3, dx3, 1, dx, 1);
    cblas_daxpy(N,  1./6, dx4, 1, dx, 1);

    free_c(dx2);
    free_c(dx3);
    free_c(dx4);
    pass_spk_id = true;

    return dx;
}


double *solve_deq_using_rk2(double* (*f) (double*, void*, void*), int N, double *x, void *arg1, void *arg2)
{
    pass_spk_id = true;
    double *dx = f(x, arg1, arg2);

    pass_spk_id = false;
    double *xtmp = (double*) malloc_c(sz_d * N);
    vdAdd(N, x, dx, xtmp);
    double *dx2 = f(xtmp, arg1, arg2);
    pass_spk_id = true;
    free_c(xtmp);

    cblas_daxpby(N, 1./2., dx2, 1, 1./2., dx, 1);
    free_c(dx2);

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
        // N++;
        
        
        // 이것만 해도 돌아감 왜? (condition: seed n*10 + 5000, n=[10, 100))
        // cells->t_fired[n][cells->num_spk[n]] = -1;
        int status_no_spike;
        double *phase = get_spike_phase(cells->num_spk[n], len, cells->t_fired[n], &status_no_spike);
        if (status_no_spike == 1){
            free(phase);
            continue;
        }

        N++;
        double *tmp_real = (double*) calloc(len, sz_d);
        double *tmp_imag = (double*) calloc(len, sz_d);
        // 왜 vdCos에 calloc이 필요한가? (malloc 에러남)

        vdCos(len, phase, tmp_real);
        vdSin(len, phase, tmp_imag); // 같이 쓰면 터짐

        vdAdd(len, tmp_real, sum_real, sum_real);
        vdAdd(len, tmp_imag, sum_imag, sum_imag);
        
        free(phase);
        free(tmp_real);
        free(tmp_imag);
    }

    if (N < 2){
        for (int i=0; i<len; i++){
            psiK[i] = 0;
            rK[i] = 0;
        }
        free(sum_real);
        free(sum_imag);
        return;
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


double *get_fft(int len, double *x_in)
{
    // x_out is the (len+2) size array
    DFTI_DESCRIPTOR_HANDLE handle = NULL;
    MKL_LONG status;

    double *x_out_tmp = (double*) malloc(sizeof(double) * (len+2));
    double *x_out = (double*) malloc(sizeof(double) * (len/2+1));

    status = DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_REAL, 1, len);
    status = DftiSetValue(handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiCommitDescriptor(handle);
    status = DftiComputeForward(handle, x_in, x_out_tmp);
    status = DftiFreeDescriptor(&handle);

    if (status && !DftiErrorClass(status, DFTI_NO_ERROR)){
        printf("Error: %s\n", DftiErrorMessage(status));
    }

    for (int n=0; n<len/2+1; n++){
        double x_real = x_out_tmp[2*n];
        double x_imag = x_out_tmp[2*n+1];
        x_out[n] = sqrt(x_real*x_real + x_imag*x_imag)/len;

        if ((n>0)  && (n<len/2)){
            x_out[n] *= 2;
        }
    }

    free(x_out_tmp);
    return x_out;
}


double *get_fft_freq(int len, double sampling_rate)
{
    double *freq = (double*) malloc(sizeof(double) * (len/2+1));
    for (int n=0; n<len/2+1; n++){
        freq[n] = sampling_rate * ((double) n) /  ((double) len);
    }

    return freq;
}


double *get_spike_phase(int n_spk, int nmax, int *nsteps, int *status_no_spike)
{
    *status_no_spike = 0;
    double *phase = (double*) calloc(nmax, sz_d);

    if (n_spk < 2){
        *status_no_spike = 1;
        return phase;
    }

    int *n0 = NULL;
    int *n1 = nsteps;
    int n_skip = t_skip_phs/_dt;
    n_spk--;

    for (int i=0; i<nmax; i++){
        if (i == *n1) {
            if (n_spk == 0){
                break;
            }
            n0 = n1;
            n1++;
            n_spk--;
        }
        
        if (n0 == NULL) continue;
        if (*n1 - *n0 > n_skip){
            phase[i] = 2*PI * (i - *n0)/((double)(*n1 - *n0));
        }
    }

    return phase;
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
        free(cells->t_fired[n]);
        cells->t_fired[n] = (int*) malloc(sz_i * _block_size);
        cells->num_spk[n] = 0;
        cells->id_fired[n] = -1;
    }
}


void append_spike(int nstep, int *num_spk, int **t_spk)
{
    (*t_spk)[*num_spk] = nstep;
    (*num_spk)++;

    // printf("append spike\n");

    if (((*num_spk) % _block_size == 0) && (*num_spk > 0)){
        // int *tmp = *t_spk;
        *t_spk = (int*) realloc(*t_spk, (*num_spk + _block_size) * sz_d);
        // printf("resize spike address %p -> %p\n", tmp, *t_spk);
        if (t_spk == NULL){
            fprintf(stderr, "There is some error in realloc while sizing %d -> %d\n", *num_spk, *num_spk+_block_size);
        }
    }
    (*t_spk)[*num_spk] = -1;
}


double *linspace(double x0, double x1, int len_x)
{
    double *x = (double*) malloc(sizeof(double) * len_x);
    if (len_x == 1){
        // printf("Too few length selected. x is set to %.2f\n", x0);
        x[0] = x0;
        return x;
    }
    for (int n=0; n<len_x; n++){
        x[n] = n*(x1-x0)/(len_x-1)+x0;
    }
    return x;
}


void downsampling(int len, double *y_org, double sample_rate, arg_t *new_vars)
{
    double fs_org = 1000./_dt;
    int nskip = fs_org/sample_rate;
    int len_new = len/nskip;

    new_vars->len = len_new;
    new_vars->x = (double*) malloc(sz_d * len_new);
    new_vars->y = (double*) malloc(sz_d * len_new);
    for (int n=0; n<len_new; n++){
        new_vars->x[n] = n*nskip*_dt;
        new_vars->y[n] = y_org[nskip*n];
    }
}


void get_fft_summary(double *vm, double sample_rate, double t_range[2], arg_t *fft_res)
{
    // t_range is the ms unit
    int n0 = t_range[0]/1000.*sample_rate;
    int len_vm = (t_range[1] - t_range[0])/1000.*sample_rate;

    double *yft_tmp = get_fft(len_vm, vm+n0);
    double *freq_tmp = get_fft_freq(len_vm, sample_rate);

    int n_cut[2];
    if (f_cut[1] > sample_rate/2) f_cut[1] = sample_rate/2;
    for (int n=0; n<2; n++) n_cut[n] = f_cut[n]/sample_rate*len_vm;

    int len_new = n_cut[1]-n_cut[0];
    fft_res->len = len_new;
    fft_res->x = (double*) malloc(sz_d * len_new);
    fft_res->y = (double*) malloc(sz_d * len_new);
    memcpy(fft_res->x, freq_tmp+n_cut[0], sz_d * len_new);
    memcpy(fft_res->y, yft_tmp+n_cut[0], sz_d * len_new);

    free(yft_tmp);
    free(freq_tmp);
}


void set_summary_rate(double fs_new){
    fs = fs_new;
}


void get_summary(int max_step, double *vm, neuron_t *cells, int *targets, res_t *summary)
{
    // downsampling
    arg_t res_vm;
    downsampling(max_step, vm, fs, &res_vm);
    summary->num_times = res_vm.len;
    summary->t = res_vm.x;
    summary->vm = res_vm.y;

    arg_t res_rk;
    double *rk_tmp, *psi_tmp;
    rk_tmp = (double*) calloc(max_step, sz_d);
    psi_tmp = (double*) calloc(max_step, sz_d);
    int flag=0;
    if (targets == NULL){
        targets = (int*) calloc(cells->num_cells, sizeof(int));
        for (int n=0; n<cells->num_cells; n++){ targets[n] = 1; }
        flag = 1;
    }
    get_Kuramoto_order_params(max_step, cells, targets, rk_tmp, psi_tmp);
    downsampling(max_step, rk_tmp, fs, &res_rk);
    summary->rk = res_rk.y;
    free(res_rk.x);
    free(rk_tmp);
    free(psi_tmp);
    if (flag == 1) free(targets);

    arg_t res_fft;
    double t_range[2] = {(max_step*_dt)-5000, (max_step-1)*_dt};
    if (t_range[0] < 0) t_range[0] = 0;
    get_fft_summary(summary->vm, fs, t_range, &res_fft);
    summary->num_freqs = res_fft.len;
    summary->freq = res_fft.x;
    summary->yf = res_fft.y;
}



void free_summary(res_t *summary)
{
    free(summary->t);
    free(summary->vm);
    free(summary->rk);
    free(summary->freq);
    free(summary->yf);
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
    free(cells->types);
}


void free_syns(syn_t *syns)
{
    free_c(syns->r);
    free_c(syns->veq);
    free_c(syns->weight);
    free_c(syns->inv_tau);

    free(syns->id_pre_neuron);
    free(syns->id_post_neuron);
    free(syns->ptr_ipost);
    free(syns->ptr_vpost);

    if (syns->type & BACKGROUND){
        free(syns->p_fire);
    } else if (syns->type & DELAY){
        free(syns->delay);
        free(syns->id_exp);
    }

    if (syns->type & STD){
        free(syns->x);
        free(syns->z);
    }
}


void free_rand_stream()
{
    vslDeleteStream(&rand_stream);
}


void destroy_mkl_buffers()
{
    mkl_free_buffers();
}


void init_network(network_info_t *info, neuron_t *cells, syn_t *syns, syn_t *bck_syns)
{
    int *cell_types = gen_types(info->num_cells, info->cell_type_ratio);
    init_cell_vars(cells, info->num_cells, info->cell_params, cell_types);

    ntk_t ntk_syn;
    init_bi_ntk(info->num_cells, info->num_cells, &ntk_syn);
    if (info->type_ntk == MEAN_DEG){
        gen_bi_random_ntk_mean_deg(cell_types, cell_types, info->mean_degs, info->gsyns, &ntk_syn);
    } else if (info->type_ntk == PROB) {
        gen_bi_random_ntk_with_type(cell_types, cell_types, info->psyns, info->gsyns, &ntk_syn);
    } else if (info->type_ntk == FIXED_IN_DEG) {
        gen_bi_random_ntk_fixed_indeg(cell_types, cell_types, info->mean_degs, info->gsyns, &ntk_syn);
    }

    init_syn_vars(syns, info->num_cells, info->type, &ntk_syn,
                        info->syn_veq, info->syn_tau, cells->v, cells->ic);
    if (info->type & DELAY){
        for (int n=0; n<syns->num_syns; n++){
            syns->delay[n] = genrand64_normal(info->t_delay_m, info->t_delay_std)/_dt;
        }
    }

    ntk_t ntk_bck;
    int *bck_types = gen_types(info->num_bck, info->bck_type_ratio);
    init_bi_ntk(info->num_bck, info->num_cells, &ntk_bck);
    gen_bi_random_ntk_with_type(bck_types, cell_types, info->pbck, info->gbck, &ntk_bck);
    init_syn_vars(bck_syns, info->num_bck, BACKGROUND, &ntk_bck, info->bck_veq, info->bck_tau, cells->v, cells->ic);
    for (int n=0; n<info->num_bck; n++){
        int tp = bck_types[n];
        bck_syns->p_fire[n] = info->frbck[tp] * _dt / 1000.;
    }

    free(cell_types);
    free(bck_types);
    free_bi_ntk(&ntk_syn);
    free_bi_ntk(&ntk_bck);
}


int *gen_types(int num, double *ratio)
{
    int *types = (int*) malloc(sizeof(int) * num);

    int n_type=0;
    int n_type_init = 0;

    for (int n=0; n<num; n++){
        if (n-n_type_init >= (double) num*ratio[n_type]){
            n_type ++;
            n_type_init = n;
        }   
        if (n_type >= MAX_TYPE){
            printf("wrong type selectd\n");
        }
        types[n] = n_type;
    }
    return types;
}

