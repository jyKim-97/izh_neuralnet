#ifndef WRITER
#define WRITER

#include "izh.h"

typedef enum _WRITER_VAR
{
    V_ONLY = 1,
    U_ONLY = 2,
    I_ONLY = 4,
    SPK_ONLY = 8

} WRITER_VAR;


typedef struct _writer_t
{
    FILE *fv;
    FILE *fu;
    FILE *fi;
    FILE *ft_spk;
    WRITER_VAR mod;
    char tag[50];

} writer_t;

void init_writer(writer_t *fid_obj, char tag[], WRITER_VAR mod);
void write_env(writer_t *fid_obj, int  num_cells, int num_bck, double tmax, double dt, int *cell_types);
void write_result(writer_t *fid_obj, int nstep, neuron_t *cells);
// void write(int nstep);
void write_data(FILE *fp, int num_x, double *x);
void write_spike(FILE *fp, int nstep, int *t_spk);
void end_writer(writer_t *fid_obj);
FILE *open_file(char fname[], char *type);
void close_file(FILE *fid);
int is_opened(writer_t *fid_obj, int mod);

#endif
