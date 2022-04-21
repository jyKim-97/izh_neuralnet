#ifndef _UTIL
#define _UTIL

#include <sys/time.h>
#include <time.h>

#define MAX_IND_NUM 10
#define _len_progbar 50
#define CHECKPOINT(tic) clock_gettime(CLOCK_MONOTONIC, &tic) 

typedef struct _progbar_t {

    int max_step;
    int div;
    struct timeval tic;
    int id_line, tot_line;
} progbar_t;


typedef struct _index_t {
    int nstep;
    int num_id;
    int len;
    int id[MAX_IND_NUM];
    int id_max[MAX_IND_NUM];
} index_t;


void set_progressbar_line(progbar_t *bar, int id_line, int tot_line);
void init_progressbar(progbar_t *bar, int max_step);
void progressbar(progbar_t *bar, int nstep);
double get_dt(struct timeval tic, struct timeval toc);
double get_dt_clock(struct timespec tic, struct timespec toc);
void print_elapsed(struct timeval start_t);
void set_index_obj(index_t *idxer, int num_index, int max_ind[]);
void next_index(index_t *idxer);
void update_index(index_t *idxer, int nstep);
// maintaining code
void print_variable(double *x, int n_print_node);
void init_multiple_lines(progbar_t *bar);
void end_multiple_lines(progbar_t *bar);

#endif