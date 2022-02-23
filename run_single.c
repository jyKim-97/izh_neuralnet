#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "parson.h"
#include "izh.h"
#include "utils.h"
#include "writer.h"

#define STP

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
    double tau_in_stp, tau_r_stp;

    double tmax;

    double sq_amp;
    double sq_inp_t[2];
    
    // name tag of the output file
    char tag[100];

} simulinfo_t;


extern double _R;
extern double _dt;
extern double default_cell_params[4][4];
extern double default_syn_veq[4];
extern double default_syn_tau[4];


void read_single_info(char fjson[100], simulinfo_t *info);
void init_simulation(simulinfo_t *info, neuron_t *cells, syn_t *syns, syn_t *bck_syns);
void run_simulation(char fjson[100]);
void gen_cell_type(int num_cells, double cell_ratio[], int *cell_types);
void alloc1d(JSON_Array *arr, double x[_n_types]);
void alloc2d(JSON_Array *arr, double x[_n_types][_n_types]);


int main(int argc, char **argv)
{

    init_random_stream(2000);

    char fname[100] = "./simul_infos/single_ntk.json";
    if (argc == 2){
        strcpy(fname, argv[1]);
    }

    run_simulation(fname);
    printf("\nDone\n");

    return 0;
}


void run_simulation(char fjson[100])
{
    _R = 0.001;

    simulinfo_t info;
    read_single_info(fjson, &info);

    neuron_t cells;
    syn_t syns;
    syn_t bck_syns;

    init_simulation(&info, &cells, &syns, &bck_syns);
    printf("write to %s\n", info.tag);

    // run simulation
    int max_step = info.tmax/_dt;
    int ninp[2] = {info.sq_inp_t[0]/_dt, info.sq_inp_t[1]/_dt};
    int num_exc = info.num_cells * info.cell_ratio[0];

    writer_t fp_obj;
    init_writer(&fp_obj, info.tag, V_ONLY | U_ONLY | I_ONLY | SPK_ONLY);
    write_env(&fp_obj, info.num_cells, info.num_bck, info.tmax, _dt, info.cell_types);

    progbar_t bar;
    init_progressbar(&bar, max_step);

    double *ic_inp = (double*) calloc(info.num_cells, sizeof(double));
    for (int i=0; i<num_exc; i++){
        ic_inp[i] = info.sq_amp;
    }

    for (int n=0; n<max_step; n++){
        double *ic;
        if ((n>ninp[0]) && (n<ninp[1])){
            ic = ic_inp;   
        } else {
            ic = NULL;
        }
        
        #ifndef STP
        update_no_delay(n, ic, &cells, &syns, &bck_syns);
        #else
        update_no_delay_stp(n, ic, &cells, &syns, &bck_syns);
        #endif

        write_result(&fp_obj, n, &cells);

        progressbar(&bar, n);
    }

    free(info.cell_types);
    free_neurons(&cells);
    free_syns(&syns);
    free_syns(&bck_syns);
    free_rand_stream();
    destroy_mkl_buffers();
}


void read_single_info(char fjson[100], simulinfo_t *info)
{
    JSON_Value *root_value = json_parse_file(fjson);
    JSON_Object *root_obj = json_value_get_object(root_value);
    JSON_Array *arr;

    info->num_cells = json_object_get_number(root_obj, "num_cells");
    _dt = json_object_get_number(root_obj, "dt");

    arr = json_object_get_array(root_obj, "cell_ratio");
    alloc1d(arr, info->cell_ratio);

    arr = json_object_get_array(root_obj, "g_syn");
    alloc2d(arr, info->g_syn);
    arr = json_object_get_array(root_obj, "p_syn");
    alloc2d(arr, info->p_syn);

    info->tau_in_stp = json_object_get_number(root_obj, "tau_in_stp");
    info->tau_r_stp  = json_object_get_number(root_obj, "tau_r_stp");

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

    #ifndef STP
    sprintf(info->tag, "%s", json_object_get_string(root_obj, "tag"));
    #else
    sprintf(info->tag, "%s_stp", json_object_get_string(root_obj, "tag"));
    #endif
    
    json_value_free(root_value);
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


void init_simulation(simulinfo_t *info, neuron_t *cells, syn_t *syns, syn_t *bck_syns)
{

    int *cell_types = (int*) malloc(sizeof(int) * info->num_cells);
    gen_cell_type(info->num_cells, info->cell_ratio, cell_types);

    init_cell_vars(cells, info->num_cells, default_cell_params, cell_types);
    info->cell_types = cell_types;

    ntk_t ntk_syn;
    init_bi_ntk(info->num_cells, info->num_cells, &ntk_syn);
    gen_bi_random_ntk_with_type(cell_types, cell_types, info->p_syn, info->g_syn, &ntk_syn);
    init_syn_vars(syns, info->num_cells, NO_DELAY, &ntk_syn, default_syn_veq,
                    default_syn_tau, cells->v, cells->ic);
    free_bi_ntk(&ntk_syn);

    syns->tau_in = info->tau_in_stp;
    syns->tau_r = info->tau_r_stp;

    ntk_t ntk_bck_syn;
    int *bck_types = (int*) calloc(info->num_bck, sizeof(int));
    double bck_veq[] = {0, 0, 0, 0};
    double bck_tau[] = {5, 5, 5, 5};

    init_bi_ntk(info->num_bck, info->num_cells, &ntk_bck_syn);
    gen_bi_random_ntk_with_type(bck_types, cell_types, &info->p_bck, &info->g_bck, &ntk_bck_syn);
    init_syn_vars(bck_syns, info->num_bck, BACKGROUND, &ntk_bck_syn, bck_veq, bck_tau,
                    cells->v, cells->ic);
    
    // set firing probability of the background synapse
    for (int n=0; n<info->num_bck; n++){
        bck_syns->p_fire[n] = info->p_fire_bck;
    }
    
    free_bi_ntk(&ntk_bck_syn);
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