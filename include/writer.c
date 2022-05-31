#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "writer.h"


double fs_write = 2000;
extern double _dt;


void init_writer(writer_t *fid_obj, char tag[], WRITER_VAR mod)
{
    char fname[100];

    fid_obj->mod = mod;
    strcpy(fid_obj->tag, tag);

    if (mod & 1){
        sprintf(fname, "%s_fv.dat", tag);
        fid_obj->fv = open_file(fname, "wb");
    }   
    if (mod & 2){
        sprintf(fname, "%s_fu.dat", tag);
        fid_obj->fu = open_file(fname, "wb");
    }
    if (mod & 4){
        sprintf(fname, "%s_fi.dat", tag);
        fid_obj->fi = open_file(fname, "wb");
    }
    if (mod & 8){
        sprintf(fname, "%s_ft_spk.txt", tag);
        fid_obj->ft_spk = open_file(fname, "w");
    }

    if (fs_write == -1){
        fid_obj->nskip = 1;
    } else {
        fid_obj->nskip = 1000./_dt/fs_write;
        if (fid_obj->nskip==0) fid_obj->nskip = 1;
    }
}


void write_summary(char tag[], res_t *simul_result)
{
    char fname[150];
    sprintf(fname, "%s_summary.info", tag);
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "num_times=%d,num_freqs=%d", simul_result->num_times, simul_result->num_freqs);
    fclose(fp);

    sprintf(fname, "%s_summary.dat", tag);
    FILE *fp2 = fopen(fname, "wb");
    fwrite(simul_result->t, sizeof(double), simul_result->num_times, fp2);
    fwrite(simul_result->vm, sizeof(double), simul_result->num_times, fp2);
    fwrite(simul_result->rk, sizeof(double), simul_result->num_times, fp2);
    fwrite(simul_result->freq, sizeof(double), simul_result->num_freqs, fp2);
    fwrite(simul_result->yf, sizeof(double), simul_result->num_freqs, fp2);
    fclose(fp2);
}


void write_result(writer_t *fid_obj, int nstep, neuron_t *cells)
{
    int num_cells = cells->num_cells;
    if (cells->nstep % fid_obj->nskip == 0){
        if (fid_obj->mod & 1){
            write_data(fid_obj->fv, num_cells, cells->v);
        }
        if (fid_obj->mod & 2){
            write_data(fid_obj->fu, num_cells, cells->u);
        }
        if (fid_obj->mod & 4){
            write_data(fid_obj->fi, num_cells, cells->ic);
        }
    }

    if (fid_obj->mod & 8){
        write_spike(fid_obj->ft_spk, nstep, cells->id_fired);
    }
}


int is_opened(writer_t *fid_obj, int mod)
{
    if (fid_obj->mod & mod){
        return 1;
    } else {
        return 0;
    }
}


FILE *open_file(char fname[], char *type)
{
    FILE *fid = fopen(fname, type);
    if (fid == NULL){
        // printf();
        char err_msg[150];
        sprintf(err_msg, "File %s is not openned", fname);
        perror(err_msg);
    }

    return fid;
}


void close_file(FILE *fid)
{
    int res = fclose(fid);
    if (res == 1){
        fprintf(stderr, "File not closed\n");
    }
}


void write_data(FILE *fp, int num_x, double *x)
{
    fwrite(x, sizeof(double), num_x, fp);
}


void write_spike(FILE *fp, int nstep, int *t_spk)
{
    int *ptr_t = t_spk;
    while (*ptr_t > -1){
        fprintf(fp, "%d-%d,", nstep, *ptr_t++);        
    }
}


void write_spike_dat(writer_t *fid_obj, neuron_t *cells)
{
    int N = cells->num_cells;

    char fname[100];
    sprintf(fname, "%s_ft_spk.info", fid_obj->tag);
    FILE *fp = fopen(fname, "w");

    int len=0;
    for (int n=0; n<N; n++){
        len += cells->num_spk[n];
        fprintf(fp, "%d,", cells->num_spk[n]);
    }
    fclose(fp);

    int id=0;
    int *t_spk_flat = (int*) malloc(sizeof(int) * len);
    for (int n=0; n<N; n++){
        int num = cells->num_spk[n];
        int *t_fired = cells->t_fired[n];
        for (int i=0; i<num; i++){
            t_spk_flat[id++] = t_fired[i];
        }
    }

    sprintf(fname, "%s_ft_spk.dat", fid_obj->tag);
    fp = fopen(fname, "wb");
    fwrite(t_spk_flat, sizeof(int), len, fp);
    fclose(fp);

    free(t_spk_flat);
}


void write_env(char tag[], network_info_t *info, int num_add, ...){
    char fname[200];
    sprintf(fname, "%s_env.json", tag);

    JSON_Value *root_value;
    JSON_Object *root_obj;

    root_value = json_value_init_object();
    root_obj = json_value_get_object(root_value);

    json_object_set_number(root_obj, "num_cells", info->num_cells);
    write_array_d(root_obj, "cell_type_ratio", info->cell_type_ratio, 2);
    write_array_i(root_obj, "mean_degs, e->", info->mean_degs[0], 2);
    write_array_i(root_obj, "mean_degs, i->", info->mean_degs[1], 2);
    write_array_d(root_obj, "g_syn, e->", info->gsyns[0], 2);
    write_array_d(root_obj, "g_syn, i->", info->gsyns[1], 2);

    write_array_d(root_obj, "g_syn, i->", info->gsyns[1], 2);
    write_array_d(root_obj, "g_syn, i->", info->gsyns[1], 2);
    
    json_object_set_number(root_obj, "num_bck", info->num_bck);
    json_object_set_number(root_obj, "fr_bck", info->frbck[0]);
    write_array_d(root_obj, "p_bck", info->pbck[0], 2);
    write_array_d(root_obj, "g_bck", info->gbck[0], 2);

    json_object_set_number(root_obj, "dt", _dt);
    json_object_set_number(root_obj, "tmax", info->tmax);
    json_object_set_number(root_obj, "fs", info->fs);
    json_object_set_number(root_obj, "type", info->type);
    json_object_set_number(root_obj, "t_delay_m", info->t_delay_m);
    json_object_set_number(root_obj, "t_delay_s", info->t_delay_std);
    json_object_set_number(root_obj, "bck_network_type", info->type_ntk);

    // read additional parameters
    va_list ap;
    va_start(ap, num_add);
    for (int n=num_add; n>0; n--)
    {
        char *var_name = va_arg(ap, char*);
        double var = va_arg(ap, double);
        json_object_set_number(root_obj, var_name, var);
    }
    va_end(ap);

    json_serialize_to_file_pretty(root_value, fname);
    json_value_free(root_value);
}

/***
void write_env(writer_t *fid_obj, int num_cells, int num_bck, double tmax, double dt, int *cell_types)
{
    char fname[200];
    sprintf(fname, "%s_env.txt", fid_obj->tag);
    FILE *fp = fopen(fname, "w");

    fprintf(fp, "num_cells=%d\n", num_cells);
    fprintf(fp, "num_bck=%d\n", num_bck);
    fprintf(fp, "tmax=%f\n", tmax);
    fprintf(fp, "dt=%f\n", dt);
    for (int n=0; n<num_cells; n++){
        fprintf(fp, "%d,", cell_types[n]);
    }
    // file list
    fprintf(fp, "\n");
    for (int i=1; i<16; i*=2){
        fprintf(fp, "%d", is_opened(fid_obj, i));
    }

    fclose(fp);
}
***/


void write_ntk(char fname[], syn_t *syns)
{
    FILE *fp = fopen(fname, "w");
    for (int n=0; n<syns->num_syns; n++){
        fprintf(fp, "%d,%d,%d,%f\n", n, syns->id_pre_neuron[n], syns->id_post_neuron[n], syns->weight[n]);
    }
    fclose(fp);
}


void write_array_d(JSON_Object *root_obj, char arr_name[], double arr1d[], int narr)
{
    json_object_set_value(root_obj, arr_name, json_value_init_array());
    JSON_Array *arr = json_object_get_array(root_obj, arr_name);
    for (int n=0; n<narr; n++){
        json_array_append_number(arr, arr1d[n]);
    }
}


void write_array_i(JSON_Object *root_obj, char arr_name[], int arr1d[], int narr)
{
    json_object_set_value(root_obj, arr_name, json_value_init_array());
    JSON_Array *arr = json_object_get_array(root_obj, arr_name);
    for (int n=0; n<narr; n++){
        json_array_append_number(arr, arr1d[n]);
    }
}


void end_writer(writer_t *fid_obj, neuron_t *cells)
{
    int mod = fid_obj->mod;
    
    if (mod & V_ONLY){
        close_file(fid_obj->fv);
    }   
    if (mod & U_ONLY){
        close_file(fid_obj->fu);
    }
    if (mod & I_ONLY){
        close_file(fid_obj->fi);
    }
    if (mod & SPK_ONLY){
        close_file(fid_obj->ft_spk);
    }
    if (mod & SPK_DAT){
        write_spike_dat(fid_obj, cells);
    }
}


