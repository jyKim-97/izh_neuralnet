#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "parson.h"
#include "izh.h"
#include "mt64.h"
#include "utils.h"

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

    double tmax;

    double sq_amp;
    double sq_inp_t[2];
    
    // name tag of the output file
    char tag[100];

} simulinfo_t;


void init_simulation(simulinfo_t *info, neuron_t *cells, syn_t *syns, bcksyn_t *bck_syns);
void gen_cell_type(int num_cells, double cell_ratio[], int *cell_types);
void read_single_info(char fjson[100], simulinfo_t *info);
void alloc1d(JSON_Array *arr, double x[_n_types]);
void alloc2d(JSON_Array *arr, double x[_n_types][_n_types]);
void open_dat_with_tag(FILE **fid, char tag[], char typename[]);
void open_file_with_tag(FILE **fid, char tag[], char typename[]);
void square_current(int n, double *ic, int num_cells, int n0, int n1, double amp);

extern double gRatio; // from Izh.c
extern double CELL_TAU[_n_types];
extern double _dt;

int main(int argc, char **argv)
{

    init_random_stream(time(NULL));

    char fname[100] = "./simul_infos/single_ntk.json";
    if (argc == 2){
        strcpy(fname, argv[1]);
    }

    // read argument
    simulinfo_t info;
    // read simulinfo params
    read_single_info(fname, &info);

    // open file writer
    neuron_t cells;
    syn_t syns;
    bcksyn_t bck_syns;

    gRatio = 0.01;
    _dt = 0.01;
    init_simulation(&info, &cells, &syns, &bck_syns);
    
    // open writer
    FILE *fv, *fu, *fi, *fr, *ft; // ft is the spike times
    open_dat_with_tag(&fv, info.tag, "v");
    open_dat_with_tag(&fu, info.tag, "u");
    open_dat_with_tag(&fi, info.tag, "i");
    open_dat_with_tag(&fr, info.tag, "r");
    open_file_with_tag(&ft, info.tag, "t");

    // write simulation info
    char fenv[100];
    sprintf(fenv, "%s_env.txt", info.tag);
    save_env(fenv, cells.num_cells, info.cell_types,
            syns.num_syns, syns.id_pre_neuron, info.tmax);

    // run simulation
    progbar_t bar;
    int max_step = info.tmax / _dt;

    int num_exc = 0;
    for (int i=0; i<info.num_cells; i++){
        if (info.cell_types[i] != 1){ num_exc++; }
    }
    free(info.cell_types);

    square_current(-1, NULL, num_exc, info.sq_inp_t[0]/_dt, info.sq_inp_t[1]/_dt, info.sq_amp);
    printf("Start single network simulation, Network Size = %d, tmax = %.1f ms, dt = %.3f ms, %d itr\n",
                info.num_cells, info.tmax, _dt, max_step);
    init_progressbar(&bar, max_step);

    write_dat(fv, cells.num_cells, cells.v);
    write_dat(fu, cells.num_cells, cells.u);
    write_dat(fi, cells.num_cells, cells.ic);
    write_dat(fr, syns.num_syns, syns.r);
    write_cell_fire(ft, cells.id_fire);
    
    for (int n=0; n<max_step; n++){
        // square input current
        square_current(n, cells.ic, 0, 0, 0, 0);
        // update all
        update_no_delay(n, &cells, &syns, &bck_syns);
        
        write_dat(fv, cells.num_cells, cells.v);
        write_dat(fu, cells.num_cells, cells.u);
        write_dat(fi, cells.num_cells, cells.ic);
        write_dat(fr, syns.num_syns, syns.r);
        write_cell_fire(ft, cells.id_fire);

        memset(cells.ic, 0, cells.num_cells*sizeof(double));

        progressbar(&bar, n);
    }
    fprintf(stderr, "\n");
    
    free_cells(&cells);
    free_syns(&syns);
    free_bcksyns(&bck_syns);

    fclose(fv);
    fclose(fu);
    fclose(fi);
    fclose(fr);
    fclose(ft);
}


void square_current(int n, double *ic, int num_cells, int n0, int n1, double amp)
{
    static int n_inp[2] = {0};
    static double amp_inp;
    static int N = 0;

    if (n == -1){
        n_inp[0] = n0;
        n_inp[1] = n1;
        amp_inp = amp;
        N = num_cells;

        return;
    }

    if ((n > n_inp[0]) && (n < n_inp[1])){
        for (int i=0; i<N; i++){
            ic[i] = amp_inp;
        }
    }
}


void open_file_with_tag(FILE **fid, char tag[], char typename[])
{
    char fname[100];
    sprintf(fname, "%s_%s.txt", tag, typename);
    *fid = fopen(fname, "w");
}


void open_dat_with_tag(FILE **fid, char tag[], char typename[])
{
    char fname[100];
    sprintf(fname, "%s_%s.dat", tag, typename);
    *fid = fopen(fname, "wb");
}


void read_single_info(char fjson[100], simulinfo_t *info)
{
    JSON_Value *root_value = json_parse_file(fjson);
    JSON_Object *root_obj = json_value_get_object(root_value);
    JSON_Array *arr;

    info->num_cells = json_object_get_number(root_obj, "num_cells");

    arr = json_object_get_array(root_obj, "cell_ratio");
    alloc1d(arr, info->cell_ratio);

    arr = json_object_get_array(root_obj, "g_syn");
    alloc2d(arr, info->g_syn);
    arr = json_object_get_array(root_obj, "p_syn");
    alloc2d(arr, info->p_syn);

    info->num_bck = json_object_get_number(root_obj, "num_bck");
    info->p_fire_bck = json_object_get_number(root_obj, "fr_bck");
    info->p_fire_bck = info->p_fire_bck / 1000 * _dt;
    info->tau_bck = json_object_get_number(root_obj, "tau_bck");

    arr = json_object_get_array(root_obj, "g_bck");
    alloc1d(arr, info->g_bck);
    arr = json_object_get_array(root_obj, "p_bck");
    alloc1d(arr, info->p_bck);

    info->tmax = json_object_get_number(root_obj, "max_time");

    info->sq_amp = json_object_get_number(root_obj, "sq_amp");
    arr = json_object_get_array(root_obj, "sq_inp_t");
    alloc1d(arr, info->sq_inp_t);

    sprintf(info->tag, "%s", json_object_get_string(root_obj, "tag"));
}


void alloc1d(JSON_Array *arr, double x[_n_types])
{
    for (int i=0; i<json_array_get_count(arr); i++){
        x[i] = json_array_get_number(arr, i);
    }
}


void alloc2d(JSON_Array *arr, double x[_n_types][_n_types])
{
    JSON_Array *arr2;

    for (int i=0; i<json_array_get_count(arr); i++){
        arr2 = json_array_get_array(arr, i);

        for (int j=0; j<json_array_get_count(arr2); j++){
            x[i][j] = json_array_get_number(arr2, j);
        }
    }   
}


void init_simulation(simulinfo_t *info, neuron_t *cells, syn_t *syns, bcksyn_t *bck_syns)
{
    ntk_t ntk_syns, ntk_bck_syns;
    
    int is_delay_on = 0;
    double *delay = NULL;

    info->cell_types = (int*) malloc(info->num_cells * sizeof(int));
    gen_cell_type(info->num_cells, info->cell_ratio, info->cell_types);

    // init cell structure
    init_cell_vars(cells, info->num_cells, info->cell_types);

    // init syn structures
    init_bi_ntk(info->num_cells, info->num_cells, &ntk_syns); // synaptic network
    gen_bi_random_ntk_with_type(info->cell_types, info->cell_types, info->p_syn, info->g_syn, &ntk_syns);
    init_syn_vars(syns, info->num_cells, is_delay_on, delay, cells->v, cells->ic, &ntk_syns);
    free_bi_ntk(&ntk_syns);

    int *bck_types = (int*) calloc(info->num_bck, sizeof(int));
    init_bi_ntk(info->num_bck, info->num_cells, &ntk_bck_syns); // background synapse
    gen_bi_random_ntk_with_type(bck_types, info->cell_types, &(info->p_bck), &(info->g_bck), &ntk_bck_syns);
    init_bcksyn_vars(bck_syns, info->num_cells, info->num_bck, info->tau_bck, &ntk_bck_syns);
    free_bi_ntk(&ntk_bck_syns);

    // set firing probability of the background synapse
    for (int n=0; n<info->num_bck; n++){
        bck_syns->fspk[n] = info->p_fire_bck;
    }
    
    free(bck_types);
}


void gen_cell_type(int num_cells, double cell_ratio[], int *cell_types)
{
    int ctp = 0;
    int num_last = num_cells*cell_ratio[0]-1;

    for (int n=0; n<num_cells; n++){
        while (n > num_last)
        {
            num_last += cell_ratio[++ctp] * num_cells;
        }
        cell_types[n] = ctp;
    }
}