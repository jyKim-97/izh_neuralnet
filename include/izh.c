#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mt64.h"
#include "mkl.h"
#include "mkl_vsl.h"
#include "izh.h"

#define PI 3.14159265359
VSLStreamStatePtr rand_stream;

const double CELL_CONST[_n_types][4] = {
            {0.02, 0.2, -65, 8},    // RS
            {0.1, 0.2, -65, 2},   // FS
            {0.02, 0.2, -55, 4},    // IB
            {0.02, 0.2, -50, 2}};   // CH

double CELL_VEQ[_n_types] = {0, -70, 0, 0};
double CELL_TAU[_n_types] = {5, 6, 5, 5};

double gRatio = 1;
double _dt = 0.005;


void init_random_stream(long int seed)
{
    vslNewStream(&rand_stream, VSL_BRNG_MT19937, seed);
    init_genrand64(seed);
}


void init_cell_vars(neuron_t *cells, int num_cells, int *cell_types)
{
    int n, ctp;
    cells->num_cells = num_cells;

    cells->v  = (double*) malloc(num_cells * sizeof(double));
    cells->u  = (double*) calloc(num_cells,  sizeof(double));
    cells->ic = (double*) calloc(num_cells,  sizeof(double));
    cells->a  = (double*) malloc(num_cells * sizeof(double));
    cells->b  = (double*) malloc(num_cells * sizeof(double));
    cells->c  = (double*) malloc(num_cells * sizeof(double));
    cells->d  = (double*) malloc(num_cells * sizeof(double));
    
    // allocate cell params
    for (n=0; n<cells->num_cells; n++)
    {
        // cells->v[n] = genrand64_real2() * 10 - 70;
        cells->v[n] = -70;
        // set parameters
        ctp = cell_types[n];
        cells->a[n] = CELL_CONST[ctp][0];
        cells->b[n] = CELL_CONST[ctp][1];
        cells->c[n] = CELL_CONST[ctp][2];
        cells->d[n] = CELL_CONST[ctp][3];
    }

    // init spike params
    cells->num_spk = (int*) calloc(num_cells,  sizeof(int));
    cells->t_spk = (int**) malloc(num_cells * sizeof(int*));
    cells->id_fire = (int*) malloc(num_cells * sizeof(int));
    for (n=0; n<num_cells; n++){
        // set t_spk & id_fire
        cells->t_spk[n]    = (int*) malloc(_block_size * sizeof(int));
        cells->t_spk[n][0] = -1;
        cells->id_fire[n]  = -1;
    }
}


void init_syn_vars(syn_t *syns, int num_cells, int is_delay_on, double *delay, double *vpost, double *ipost, ntk_t *ntk)
{
    // set type of the synapse
    syns->is_delay_on = is_delay_on;
    
    syns->num_cells = num_cells;
    syns->num_pre = num_cells;

    // get the # of syns
    int n, i, num_syns=0;
    for (n=0; n<ntk->num_pre; n++){
        num_syns += ntk->num_edges[n];
    }
    syns->num_syns = num_syns;

    // create objs
    syns->R = gRatio;
    syns->r           = (double*) calloc(num_syns, sizeof(double));
    syns->weight      = (double*) malloc(num_syns * sizeof(double));
    syns->veq         = (double*) malloc(num_syns * sizeof(double));
    syns->inv_tau     = (double*) malloc(num_syns * sizeof(double));
    syns->ptr_vpost   = (double**) malloc(num_syns * sizeof(double*));
    syns->ptr_ipost   = (double**) malloc(num_syns * sizeof(double*));
    syns->id_pre_neuron   = (int*) malloc(num_syns * sizeof(int));
    syns->id_post_neuron   = (int*) malloc(num_syns * sizeof(int));

    if (is_delay_on == 1){
        syns->delay    = (int*) calloc(num_syns, sizeof(int));
        syns->id_pre   = (int*) malloc(num_syns * sizeof(int));
        syns->id_t_spk = (int*) calloc(num_syns, sizeof(int));
    } else { // used to update firing
        syns->n_pre2syn  = (int*) calloc(syns->num_pre,  sizeof(int));
        syns->id_pre2syn = (int**) malloc(syns->num_pre * sizeof(int*));
    }

    // set variables
    int pre_type, id_syn=0, n_post;

    for (n=0; n<ntk->num_pre; n++){
        pre_type = ntk->node_types[n];
        for (i=0; i<ntk->num_edges[n]; i++){
            syns->veq[id_syn] = CELL_VEQ[pre_type];
            syns->inv_tau[id_syn] = 1/CELL_TAU[pre_type];
            syns->id_pre_neuron[id_syn] = n; 

            n_post = ntk->adj_list[n][i];
            syns->ptr_vpost[id_syn] = vpost + n_post;
            syns->ptr_ipost[id_syn] = ipost + n_post;
            syns->id_post_neuron[id_syn] = n_post;

            id_syn++;
        }
    }

    // set id
    if (is_delay_on == 0){
        map_pre2syn(syns, ntk);
    }

    // allocate synapse weight
    id_syn=0;
    for (n=0; n<num_cells; n++){
        for (i=0; i<ntk->num_edges[n]; i++) {
            syns->weight[id_syn++] = ntk->strength[n][i];
        }
    }
}


void init_bcksyn_vars(bcksyn_t *bck_syns, int num_cells, int num_bck, double tau_bck, ntk_t *ntk)
{
    bck_syns->num_cells = num_cells;
    bck_syns->num_bck   = num_bck;

    bck_syns->R = gRatio;
    bck_syns->r       = (double*) calloc(num_bck, sizeof(double));
    bck_syns->weight  = (double*) calloc(num_bck * num_cells, sizeof(double));
    bck_syns->inv_tau = (double*) malloc(num_bck * sizeof(double));
    bck_syns->fspk    = (double*) malloc(num_bck * sizeof(double));
    bck_syns->syn_act = (double*) calloc(num_bck, sizeof(double));

    // set variables
    int n;
    for (n=0; n<num_bck; n++){
        bck_syns->inv_tau[n] = 1/tau_bck;
    }

    // allocate synapse weight
    int i,id_pre,id_post;
    for (id_pre=0; id_pre<num_bck; id_pre++){
        for (i=0; i<ntk->num_edges[id_pre]; i++){
            id_post = ntk->adj_list[id_pre][i];
            // column-major matrix, row is the post
            n = id_pre * num_cells + id_post;
            bck_syns->weight[n] = ntk->strength[id_pre][i];
        }
    }
}


void update_no_delay(int nstep, neuron_t *cells, syn_t *syns, bcksyn_t *bck_syns)
{   
    #ifndef _NO_EXT_SYN
    update_bcksyns(bck_syns);
    add_bcksyn_current(bck_syns, cells);
    #endif

    add_syn_current(syns, cells);
    update_cells(nstep, cells);
    update_syns_no_delay(nstep, syns, cells->id_fire);
}


void update_delay(int nstep, neuron_t *cells, syn_t *syns, bcksyn_t *bck_syns)
{
    #ifndef _NO_EXT_SYN
    update_bcksyns(bck_syns);
    add_bcksyn_current(bck_syns, cells);
    #endif

    add_syn_current(syns, cells);
    update_cells(nstep, cells);
    update_syns_delay(nstep, syns, cells->t_spk);
}


void update_cells(int nstep, neuron_t *cells)
{
    int n, N = cells->num_cells;
    double *dv, *du;

    dv = (double*) malloc(N * sizeof(double));
    du = (double*) malloc(N * sizeof(double));
    solve_deq_using_rk4(f_dv, N, dv, cells->v, cells, NULL);
    solve_deq_using_rk4(f_du, N, du, cells->u, cells, NULL);
    cblas_daxpy(N, 1, dv, 1, cells->v, 1);
    cblas_daxpy(N, 1, du, 1, cells->u, 1);

    // check spike
    double *ptr_v = cells->v;
    int *ptr_id = cells->id_fire;

    for (n=0; n<N; n++){
        if (*ptr_v > 30){
            *ptr_id++ = n;
            *ptr_v = cells->c[n];
            cells->u[n] += cells->d[n];
        }
        ptr_v++;
    }
    *ptr_id = -1;

    // append spike to t_spk
    ptr_id = cells->id_fire;
    while (*ptr_id > -1){
        append_spike(nstep, cells->num_spk + *ptr_id, cells->t_spk + *ptr_id);
        ptr_id++;
    }

    free(dv); free(du);
}


void update_syns_no_delay(int nstep, syn_t *syns, void *id_fire_tmp)
{
    int *id_fire = (int*) id_fire_tmp;
    int num_syn=syns->num_syns;
    double *dr = (double*) malloc(num_syn * sizeof(double));

    solve_deq_using_rk4(f_dr_no_delay, num_syn, dr, syns->r, syns, id_fire);

    cblas_daxpy(num_syn, 1, dr, 1, syns->r, 1);

    free(dr);
}


void update_syns_delay(int nstep, syn_t *syns, int **t_spk)
{   
    int num_syns=syns->num_syns, n;

    // get activated synapse
    double *syn_act = (double*) calloc(num_syns, sizeof(double));
    double R = syns->R;
    int t0, dt;

    // #pragma omp parallel for
    for (n=0; n<num_syns; n++){
        int *ptr_spk = syns->id_t_spk + n;

        t0 = t_spk[syns->id_pre[n]][*ptr_spk];
        dt = nstep - t0 - syns->delay[n];

        if (dt == 0){
            syn_act[n] = R;
            (*ptr_spk)++;
        }
    }

    // update r
    double *dr = (double*) malloc(num_syns * sizeof(double));

    solve_deq_using_rk4(f_dr_delay, syns->num_syns, dr, syns->r, syns, syn_act);
    cblas_daxpy(syns->num_syns, 1., dr, 1, syns->r, 1);

    free(syn_act);
    free(dr);
}


void update_bcksyns(bcksyn_t *bck_syns)
{

    int n, N=bck_syns->num_bck;
    // gen spike
    double *parr = (double*) malloc(sizeof(double)*N);
    double *ptr_p;
    // double p;

    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rand_stream, N, parr, 0, 1);
    ptr_p = parr;

    memset(bck_syns->syn_act, 0, N*sizeof(double));
    for (n=0; n<N; n++){
        // p = genrand64_real2();
        if (*ptr_p++ < bck_syns->fspk[n]){
            bck_syns->syn_act[n] = 1;
        }
    }

    free(parr);

    // update state variable - r
    double *dr = (double*) malloc(N * sizeof(double));

    solve_deq_using_rk4(f_dr_bck, N, dr, bck_syns->r, bck_syns, bck_syns->syn_act);
    cblas_daxpy(N, 1., dr, 1, bck_syns->r, 1);

    free(dr);
}


void add_bcksyn_current(bcksyn_t *bck_syns, neuron_t *cells)
{
    int num_cells=bck_syns->num_cells;
    double *isyn = (double*) calloc(num_cells, sizeof(double));

    cblas_dgemv(CblasColMajor, CblasNoTrans, num_cells, bck_syns->num_bck,
                1., bck_syns->weight, num_cells, bck_syns->r, 1, 0, isyn, 1);

    vdMul(num_cells, cells->v, isyn, isyn);
    cblas_daxpy(num_cells, -1, isyn, 1, cells->ic, 1);
    free(isyn);
}


void add_syn_current(syn_t *syns, neuron_t *cells)
{
    int num_syns=syns->num_syns, n;

    double *isyn = (double*) malloc(num_syns * sizeof(double));
    
    // copy vpost to isyn
    for (n=0; n<num_syns; n++){
        isyn[n] = *(syns->ptr_vpost[n]);
    }

    // Isyn = weight .* r .* (vpost - veq)
    cblas_daxpy(num_syns, -1, syns->veq, 1, isyn, 1);
    vdMul(num_syns, syns->r, isyn, isyn);
    vdMul(num_syns, syns->weight, isyn, isyn);

    // add current to each post-synaptic neuron
    for (n=0; n<num_syns; n++){
        *(syns->ptr_ipost[n]) -= isyn[n];
    }

    free(isyn);
}


void solve_deq_using_rk4(void (*f) (double*, double*, void*, void*), int N, double *dx, double *x, void *arg1, void *arg2)
{
    double *dx2, *dx3, *dx4, *xnew;

    int sz = sizeof(double) * N;

    // init array, sz is the length of the array
    xnew = (double*) malloc(sz);
    dx2  = (double*) malloc(sz);
    dx3  = (double*) malloc(sz);
    dx4  = (double*) malloc(sz);

    // solve equation using RK4 mehtod
    f(dx, x, arg1, arg2); // 1st step, dx = dx1

    memcpy(xnew, x, sz);   // 2nd step
    cblas_daxpy(N, 0.5, dx, 1, xnew, 1);
    f(dx2, x, arg1, arg2);

    memcpy(xnew, x, sz);   // 3rd step
    cblas_daxpy(N, 0.5, dx2, 1, xnew, 1);
    f(dx3, x, arg1, arg2);

    memcpy(xnew, x, sz);   // 4th step
    cblas_daxpy(N, 1., dx3, 1, xnew, 1);
    f(dx4, x, arg1, arg2);

    // add k1/6 + k2/3 +  k3/3 + k4/6
    cblas_daxpy(N, 1/3, dx2, 1/6, dx, 1);
    cblas_daxpy(N, 1/3, dx3, 1,   dx, 1);
    cblas_daxpy(N, 1/6, dx3, 1,   dx, 1);

    // end objs
    free(xnew);
    free(dx2); free(dx3); free(dx4);
}


void f_dv(double *dv, double *v, void *arg_cells, void *arg_null)
{
    neuron_t *cells = (neuron_t*) arg_cells;
    int N = cells->num_cells;

    vdMul(N, v, v, dv);
    cblas_dscal(N, 0.04, dv, 1);
    cblas_daxpy(N, 5., v, 1, dv, 1);
    cblas_daxpy(N, -1., cells->u, 1, dv, 1);
    cblas_daxpy(N, 1., cells->ic, 1, dv, 1);

    for (int n=0; n<cells->num_cells; n++){ dv[n] += 140; }

    cblas_dscal(N, _dt, dv, 1);
}


void f_du(double *du, double *u, void *arg_cells, void *arg_null)
{
    neuron_t *cells = (neuron_t*) arg_cells;
    int N = cells->num_cells;
    
    vdMul(N, cells->b, cells->v, du); // du = b.*v
    cblas_daxpy(N, -1., u, 1, du, 1); // du = (b.*v - u)
    vdMul(N, cells->a, du, du);       // du = a.*(b.*v - u)
    cblas_dscal(N, _dt, du, 1);
}


void f_dr_delay(double *dr, double *r, void *arg_syns, void *arg_syn_act)
{
    syn_t *syns = (syn_t*) arg_syns;
    double *syn_act = (double*) arg_syn_act;

    int num_syns = syns->num_syns;

    memcpy(dr, r, num_syns * sizeof(double));
    // cblas_dcopy(num_syns, r, 1, dr, 1);
    cblas_dscal(num_syns, -1., dr, 1);
    cblas_daxpy(num_syns, syns->R, syn_act, 1, dr, 1);
    vdMul(num_syns, dr, syns->inv_tau, dr);
    cblas_dscal(num_syns, _dt, dr, 1);
}


void f_dr_no_delay(double *dr, double *r, void *arg_syns, void *arg_id_fire)
{
    syn_t *syns = (syn_t*) arg_syns;
    int *id_fire = (int*) arg_id_fire;
    int i, id_syn, *m = id_fire, num_syns=syns->num_syns;

    // cblas_dcopy(num_syns, r, 1, dr, 1);
    // memset(dr, 0, sizeof(double) * num_syns);
    memcpy(dr, r, num_syns * sizeof(double));
    cblas_dscal(num_syns, -1., dr, 1);

    while (*m != -1){
        for (i=0; i<syns->n_pre2syn[*m]; i++){
            id_syn = syns->id_pre2syn[*m][i];
            dr[id_syn] += syns->R;
        }
        m++;
    }

    vdMul(num_syns, syns->inv_tau, dr, dr);
    cblas_dscal(num_syns, _dt, dr, 1);
}


void f_dr_bck(double *dr, double *r, void *arg_syns, void *arg_syn_act)
{
    bcksyn_t *bck_syns = (bcksyn_t*) arg_syns;
    double *syn_act = (double*) arg_syn_act;
    int N = bck_syns->num_bck;

    memcpy(dr, r, N*sizeof(double));
    // cblas_dcopy(N, r, 1, dr, 1);    
    cblas_dscal(N, -1., dr, 1);
    cblas_daxpy(N, bck_syns->R, syn_act, 1, dr, 1);
    vdMul(N, bck_syns->inv_tau, dr, dr);
    cblas_dscal(N, _dt, dr, 1);
}


void map_pre2syn(syn_t *syns, ntk_t *ntk)
{
    int id_syn=0, n, i;

    if (syns->is_delay_on == 1){
        for (n=0; n<ntk->num_pre; n++){
            for (i=0; i<ntk->num_edges[n]; i++){
                syns->id_pre[id_syn++] = n;
            }   
        }
    } else { // delay == 0
        memcpy(syns->n_pre2syn, ntk->num_edges, syns->num_pre * sizeof(int));

        for (n=0; n<ntk->num_pre; n++){
            (syns->id_pre2syn)[n] = (int*) malloc(ntk->num_edges[n] * sizeof(int));
            for (i=0; i<ntk->num_edges[n]; i++){
                syns->id_pre2syn[n][i] = id_syn++;
            }
        }   
    }
}

// get Kuramoto order parameters
void get_Kuramoto_order_params(int len, neuron_t *cells, double *rK, double *psiK)
{   
    int n, i, num_cells=cells->num_cells;
    double *phase, *sum_real, *sum_imag, *ptr_sr, *ptr_si, *ptr_p;

    phase    = (double*) malloc(sizeof(double) * len);
    sum_real = (double*) malloc(sizeof(double) * len);
    sum_imag = (double*) malloc(sizeof(double) * len);

    for (n=0; n<num_cells; n++){
        get_spike_phase(cells->num_spk[n], len, cells->t_spk[n], phase);

        ptr_sr = sum_real;
        ptr_si = sum_imag;
        ptr_p  = phase;

        for (i=0; i<len; i++){
            *ptr_sr++ += cos(*ptr_p);
            *ptr_si++ += sin(*ptr_p++);
        }
    }

    cblas_dscal(len, 1./num_cells, sum_real, 1);
    cblas_dscal(len, 1./num_cells, sum_imag, 1);

    // get Kuramoto order parameter
    ptr_sr = sum_real;
    ptr_si = sum_imag;
    ptr_p  = psiK;
    for (i=0; i<len; i++){
        *ptr_p++ = atan2(*ptr_si++, *ptr_sr++); // +- phi
    }

    cblas_dnrm2(len, sum_real, 1);

    ptr_sr = sum_real;
    ptr_si = sum_imag;
    ptr_p  = rK;
    for (i=0; i<len; i++){
        *ptr_p++ = sqrt(pow(*ptr_sr++, 2) + pow(*ptr_si++, 2));
    }

    free(phase);
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


void append_spike(int nstep, int *id, int **t_spk)
{
    (*t_spk)[*id] = nstep;
    (*id)++;
    if (((*id) % _block_size == 0) && (*id > 0)){
        *t_spk = (int*) realloc(*t_spk, (*id + _block_size) * sizeof(double));
    }
    (*t_spk)[*id] = -1;
}


double get_avg(int num_x, double *x)
{   
    // get average of the x
    double x_avg = 0;
    for (int n=0; n<num_x; n++){
        x_avg += x[n];
    }
    x_avg /= num_x;
    
    // x_avg /= num_x;
    return x_avg;
}


void reset_spike(neuron_t *cells)
{
    for (int n=0; n<cells->num_cells; n++){
        for (int i=0; i<cells->num_spk[n]; i++){
            cells->t_spk[n][i] = 0;
        }
        cells->num_spk[n] = 0;
    }
}


void save_env(char fname[100], int num_cells, int *cell_types, int num_syns, int *id_presyns, double tmax)
{
    FILE *fid = fopen(fname, "w");
    fprintf(fid, "%d, %d, %f, %f\n", num_cells, num_syns, tmax, _dt);
    // save cell types
    int n;
    fprintf(fid, "%d", cell_types[0]);
    for (n=1; n<num_cells; n++){
        fprintf(fid, ",%d", cell_types[n]);
    }
    fprintf(fid, "\n");

    // save pre neuron ids
    fprintf(fid, "%d", id_presyns[0]);
    for (n=1; n<num_syns; n++){
        fprintf(fid, ",%d", id_presyns[n]);
    }
    fprintf(fid, "\n");

    fclose(fid);
}


void write_dat(FILE *fid, int N, double *x)
{
    fwrite(x, sizeof(double), N, fid);
}


void write_cell_fire(FILE *fid, int *id_fire)
{
    int *id = id_fire;
    while (*id != -1){
        fprintf(fid, "%d,", *(id++));
    }
    fprintf(fid, "\n");
}


void free_cells(neuron_t *cells)
{
    free(cells->v);
    free(cells->u);
    free(cells->ic);
    free(cells->a);
    free(cells->b);
    free(cells->c);
    free(cells->d);
    free(cells->num_spk);
    free(cells->id_fire);

    for (int n=0; n<cells->num_cells; n++){
        free(cells->t_spk[n]);
    }
    free(cells->t_spk);
}


void free_syns(syn_t *syns)
{
    int n;

    free(syns->r);
    free(syns->weight);
    free(syns->veq);
    free(syns->inv_tau);
    free(syns->ptr_vpost);
    free(syns->ptr_ipost);
    free(syns->id_pre_neuron);
    free(syns->id_post_neuron);

    if (syns->is_delay_on == 1){
        free(syns->delay);
        free(syns->id_pre);
        free(syns->id_t_spk);
    } else {
        free(syns->n_pre2syn);
        for (n=0; n<syns->num_pre; n++){
            free(syns->id_pre2syn[n]);
        }
        free(syns->id_pre2syn);
    }
}


void free_bcksyns(bcksyn_t *bck_syns)
{
    free(bck_syns->r);
    free(bck_syns->weight);
    free(bck_syns->inv_tau);
    free(bck_syns->fspk);
    free(bck_syns->syn_act);
}
