#ifndef WRITER
#define WRITER

#include "izh.h"
#include "parson.h"

typedef enum _WRITER_VAR
{
    V_ONLY = 1,
    U_ONLY = 2,
    I_ONLY = 4,
    SPK_ONLY = 8,
    SPK_DAT = 16

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
void write_env(char tag[], network_info_t *info, int num_add, ...);
void write_summary(char tag[], res_t *simul_result);
void write_ntk(char fname[], syn_t *syns);
// void write_env(writer_t *fid_obj, int  num_cells, int num_bck, double tmax, double dt, int *cell_types);
void write_result(writer_t *fid_obj, int nstep, neuron_t *cells);
void write_data(FILE *fp, int num_x, double *x);
void write_spike(FILE *fp, int nstep, int *t_spk);
void write_spike_dat(writer_t *fid_obj, neuron_t *cells);
FILE *open_file(char fname[], char *type);
void end_writer(writer_t *fid_obj);
void close_file(FILE *fid);
int is_opened(writer_t *fid_obj, int mod);
void write_array_d(JSON_Object *root_obj, char arr_name[], double arr1d[], int narr);
void write_array_i(JSON_Object *root_obj, char arr_name[], int arr1d[], int narr);

#endif
