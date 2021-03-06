#include <stdio.h>
#include "utils.h"
#include <stdbool.h>


static bool print_multi_line = false;


void set_progressbar_line(progbar_t *bar, int id_line, int tot_line)
{
    bar->id_line = id_line;
    bar->tot_line = tot_line;
    print_multi_line = true;
}


void init_multiple_lines(progbar_t *bar)
{
    for (int n=bar->tot_line-1; n!=0; n--){
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\r\033[%dB", bar->tot_line-1);
}


void end_multiple_lines(progbar_t *bar)
{
    for (int n=bar->tot_line; n!=0; n--){
        fprintf(stderr, "\n");
    }
}


void init_progressbar(progbar_t *bar, int max_step)
{
    bar->max_step = max_step;
    bar->div = max_step/_len_progbar;
    gettimeofday(&(bar->tic), NULL);
}


void progressbar(progbar_t *bar, int nstep)
{
    if ((nstep+1)%(bar->div) == 0){
        int percent = (double)nstep/(bar->max_step)*100;
        fprintf(stderr, "\r\033[%dB[", bar->id_line);
        int i, nbar = (nstep+1)/(bar->div);
        for (i=0; i<nbar; i++){
            fprintf(stderr, "=");
        }
        for (i=nbar; i<_len_progbar; i++){
            fprintf(stderr, " ");
        }
        fprintf(stderr, "]");
        // get time
        struct timeval toc;
        double dt, pred_end_time;

        gettimeofday(&toc, NULL);
        dt = get_dt(bar->tic, toc);
        pred_end_time = dt * ((double) bar->max_step) / ((double) nstep);

        fprintf(stderr, "(%5.1fs / %5.1fs)", dt, pred_end_time);
        fprintf(stderr, "\r[%3d%%", percent);
        fprintf(stderr, "\033[%dA", bar->id_line);
    }
}


double get_dt(struct timeval tic, struct timeval toc)
{
    int sec, msec;

    sec = toc.tv_sec - tic.tv_sec;
    msec = (toc.tv_usec - tic.tv_usec)/1e3;

    if (msec < 0){
        sec -= 1;
        msec += 1e3;
    }

    double dt = sec + ((double) msec) * 1e-3;

    return dt;
}


double get_dt_clock(struct timespec tic, struct timespec toc)
{
    long sec, nsec;

    sec = toc.tv_sec - tic.tv_sec;
    nsec = (toc.tv_nsec - tic.tv_nsec);

    double dt = sec + ((double) nsec) * 1e-9;
    return dt;
}


void print_elapsed(struct timeval start_t)
{
    int sec, msec, usec, x;
    struct timeval end_t;

    gettimeofday(&end_t, NULL);


    sec = end_t.tv_sec - start_t.tv_sec;
    usec = end_t.tv_usec - start_t.tv_usec;
    x = usec / 1e3;
    msec = x;
    usec -= x * 1e3;

    if (usec < 0){
        msec -= 1;
        usec += 1e3;
    }
    if (msec < 0){
        sec -= 1;
        msec += 1e3;
    }

    printf("elapsed time = %ds %dms %dus\n", sec, msec, usec);

}


void print_variable(double *x, int n_print_node)
{
    for (int n=0; n<n_print_node; n++){
        fprintf(stderr, "%5.2f, ", x[n]);
    }
    fprintf(stderr, "\n");
}


void set_index_obj(index_t *idxer, int num_index, int max_ind[]){
    idxer->nstep = 0;
    idxer->num_id = num_index;
    idxer->len = 1;
    for (int n=0; n<num_index; n++){
        idxer->id[n] = 0;
        idxer->id_max[n] = max_ind[n];
        idxer->len *= max_ind[n];
    }
}


void next_index(index_t *idxer){

    if (idxer->nstep+1 == idxer->len){
        printf("Index out of range\n");
        return;
    }

    int nid = idxer->num_id-1;
    idxer->nstep++;
    idxer->id[nid]++;
    while (idxer->id[nid] >= idxer->id_max[nid]){
        idxer->id[nid] -= idxer->id_max[nid];
        idxer->id[nid-1]++;
        nid--;
    }
}


void update_index(index_t *idxer, int nstep){
    if (nstep >= idxer->len){
        printf("Out of range\n");
        return;
    }

    idxer->nstep = nstep;
    int div = idxer->len;
    for (int n=0; n<idxer->num_id; n++){
        div /= idxer->id_max[n];
        idxer->id[n] = nstep / div;
        nstep -= idxer->id[n] * div;
    }
}
