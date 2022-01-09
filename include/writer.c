#include <stdio.h>
#include <string.h>
#include <writer.h>


void init_writer(writer_t *fid, char tag[], WRITER_VAR mod)
{
    char fname[200];

    if (mod & 1){
        sprintf(fname, "%s_v.dat", tag);
        fid->v = fopen(fname, "wb");
    }   
    if (mod & 2){
        sprintf(fname, "%s_u.dat", tag);
        fid->u = fopen(fname, "wb");
    }
    if (mod & 4){
        sprintf(fname, "%s_i.dat", tag);
        fid->i = fopen(fname, "wb");
    }
    if (mod & 8){
        sprintf(fname, "%s_t_spk.txt", tag);
        fid->t_spk = fopen(fname, "w");
    }
}




