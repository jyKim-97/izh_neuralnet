#ifndef _UTIL
#define _UTIL

#include <sys/time.h>
#include <time.h>

#define _len_progbar 50
#define CHECKPOINT(tic) clock_gettime(CLOCK_MONOTONIC, &tic) 

typedef struct _progbar_t {

    int max_step;
    int div;
    struct timeval tic;
    int id_line;

} progbar_t;

void set_progressbar_line(progbar_t *bar, int id_line);
void init_progressbar(progbar_t *bar, int max_step);
void progressbar(progbar_t *bar, int nstep);
double get_dt(struct timeval tic, struct timeval toc);
double get_dt_clock(struct timespec tic, struct timespec toc);
void print_elapsed(struct timeval start_t);
// maintaining code
void print_variable(double *x, int n_print_node);

#endif