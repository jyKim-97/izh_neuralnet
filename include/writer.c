#include <stdio.h>
#include <string.h>
#include "writer.h"


void init_writer(writer_t *fid_obj, char tag[], WRITER_VAR mod)
{
    char fname[200];

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
}


void write(writer_t *fid_obj, int nstep, neuron_t *cells)
{
    int num_cells = cells->num_cells;
    if (fid_obj->mod & 1){
        write_data(fid_obj->fv, num_cells, cells->v);
    }
    if (fid_obj->mod & 2){
        write_data(fid_obj->fu, num_cells, cells->u);
    }
    if (fid_obj->mod & 4){
        write_data(fid_obj->fi, num_cells, cells->ic);
    }
    if (fid_obj->mod & 8){
        write_spike(fid_obj->ft_spk, nstep, cells->id_fired);
    }
}


void write_env(writer_t *fid_obj, int num_cells, int num_bck, double tmax, double dt, int *cell_types)
{
    char fname[100];
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


int is_opened(writer_t *fid_obj, int mod)
{
    if (fid_obj->mod & mod){
        return 1;
    } else {
        return 0;
    }
}


FILE *open_file(char fname[100], char *type)
{
    FILE *fid = fopen(fname, type);
    if (fid == NULL){
        // printf();
        char err_msg[100];
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


void write_data (FILE *fp, int num_x, double *x)
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


void end_writer(writer_t *fid_obj)
{
    int mod = fid_obj->mod;
    
    if (mod & 1){
        close_file(fid_obj->fv);
    }   
    if (mod & 2){
        close_file(fid_obj->fu);
    }
    if (mod & 4){
        close_file(fid_obj->fi);
    }
    if (mod & 8){
        close_file(fid_obj->ft_spk);
    }
}


