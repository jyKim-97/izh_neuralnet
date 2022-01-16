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
        // fid_obj->fv = fopen(fname, "wb");
        // check_open(fid_obj->fv, "v");
    }   
    if (mod & 2){
        sprintf(fname, "%s_fu.dat", tag);
        fid_obj->fu = open_file(fname, "wb");
        // fid_obj->fu = fopen(fname, "wb");
        // check_open(fid_obj->fu, "u");
    }
    if (mod & 4){
        sprintf(fname, "%s_fi.dat", tag);
        fid_obj->fi = open_file(fname, "wb");
        // fid_obj->fi = fopen(fname, "wb");
        // check_open(fid_obj->fi, "ic");
    }
    if (mod & 8){
        sprintf(fname, "%s_ft_spk.txt", tag);
        fid_obj->ft_spk = open_file(fname, "w");
        // fid_obj->ft_spk = fopen(fname, "w");
        // check_open(fid_obj->ft_spk, "spike");
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


