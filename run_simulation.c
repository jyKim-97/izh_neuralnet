#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "izh.h"
#include "utils.h"
#include "writer.h"
#include "mt64.h"

extern double _dt, fs;
extern double default_syn_veq[MAX_TYPE];
extern double default_syn_tau[MAX_TYPE];


char env_tag[100] = "test";
long seed = 1000;


void run(double tmax, WRITER_VAR save_var_types);
void init_simulation(network_info_t *info, neuron_t *cells, syn_t *syns, syn_t *bck_syns);
void set_cell_params(neuron_t *cells);
void set_env(network_info_t *info);


int main(int argc, char **argv){

    if (argc == 2){
        sprintf(env_tag, "%s", argv[1]);
    }
    init_random_stream(seed);
    run(5000, V_ONLY|SPK_DAT);
    free_rand_stream();
    return 0;
}


void run(double tmax, WRITER_VAR save_var_types){
    neuron_t cells;
    syn_t syns, bck_syns;
    network_info_t info={0,};

    set_env(&info);
    init_simulation(&info, &cells, &syns, &bck_syns);
    info.tmax = tmax;

    int max_step = tmax/_dt;

    writer_t fp_obj;
    init_writer(&fp_obj, env_tag, save_var_types);
    write_env(env_tag, &info, 1, "random_seed", (double) seed);
    
    printf("Run simulation -> %s\n", env_tag);
    progbar_t bar;
    init_progressbar(&bar, max_step);

    double *vm = (double*) malloc(sizeof(double) * max_step);
    int num_exc = cells.num_cells * 0.8;

    for (int n=0; n<max_step; n++){
        update(n, NULL, &cells, &syns, &bck_syns);
        
        for (int i=0; i<num_exc; i++) vm[n] += cells.v[i];
        vm[n] /= num_exc;

        write_result(&fp_obj, n, &cells);
        progressbar(&bar, n);
    }
    if (save_var_types | SPK_DAT) write_spike_dat(&fp_obj, &cells);

    res_t summary;
    get_summary(max_step, vm, &cells, NULL, &summary);
    write_summary(env_tag, &summary);
    free_summary(&summary);

    end_writer(&fp_obj);
    free_neurons(&cells);
    free_syns(&syns);
    free_syns(&bck_syns);
    destroy_mkl_buffers();

    printf("[Done\n");
}


void init_simulation(network_info_t *info, neuron_t *cells, syn_t *syns, syn_t *bck_syns){
    init_network(info, cells, syns, bck_syns);
    set_cell_params(cells);
}


void set_cell_params(neuron_t *cells){
    // gen heterogeneous network
    for (int n=0; n<cells->num_cells; n++){
        int ctp = cells->types[n];
        if (ctp == 0){
            double sgm = genrand64_real2();
            cells->a[n] = 0.02;
            cells->b[n] = 0.2;
            cells->c[n] = -65 + 15 * sgm*sgm;
            cells->d[n] = 8 - 6 * sgm*sgm;
        } else if (ctp == 1){
            cells->a[n] = 0.02;
            cells->b[n] = 0.25;
            cells->c[n] = -65;
            cells->d[n] = 2;
        } else if (ctp > 1) {
            printf("Unexpected cell type included, id=%d, type=%d\n", n, ctp);
        }
    }
}


void set_env(network_info_t *info){

    double g_const = 0.001;
    _dt = 0.05;
    info->type_p = NONE;

    info->num_cells = 500;
    info->cell_type_ratio[0] = 0.8;
    info->cell_type_ratio[1] = 0.2;
    info->mean_degs[0][0] = 80;
    info->mean_degs[0][1] = 20;
    info->mean_degs[1][0] = 120;
    info->mean_degs[1][1] = 30;

    info->gsyns[0][0] = 20*g_const;
    info->gsyns[0][1] = 20*g_const;
    info->gsyns[1][0] = 20*g_const;
    info->gsyns[1][1] = 20*g_const;

    info->num_bck = 1000;
    info->bck_type_ratio[0] = 1;
    info->frbck[0] = 1000; // hz
    info->pbck[0][0] = 0.1;
    info->pbck[0][1] = 0.1;
    info->gbck[0][0] = 1.5*g_const;
    info->gbck[0][1] = 1.5*g_const;

    for (int i=0; i<2; i++){
        info->syn_tau[i] = default_syn_tau[i];
        info->syn_veq[i] = default_syn_veq[i];
        info->bck_tau[i] = default_syn_tau[i];
        info->bck_veq[i] = default_syn_veq[i];
    }

    info->fs = fs;
    info->t_delay_m = 0;
    info->t_delay_std = 0;
}