#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parson.h"
#include "easyIzh.h"
#include "utils.h"


/***
 * Library that help easier to use
 * use json file to pass simulation information
 * format of the .json file (tree)
 * 
 * TODO: 채우기


***/


extern double gRatio;


void read_simulinfo(simul_info_t *info, char jsonfname[100])
{
    
}


void run_simulation(simul_info_t *info)
{
    // init structures
    neuron_t cells;
    syn_t syns;
    bcksyn_t bck_syns;

    // initializing
    gRatio = 0.05;
    init_simulation(info, &cells, &syns, &bck_syns);

    // run simulation
    progbar_t bar;
    int max_step = info->tmax/_dt;
    writer_t fdat;

    printf("Start simulation, Network Size = %d, tmax = %.1f ms, dt = %.2f ms, %d itr\n",
            info->num_cells, info->tmax, _dt, max_step);
    init_writer(&fdat, "./data/testnew", info->num_cells, info->cell_types, info->tmax);
    free(info->cell_types);

    int n0 = (100. / _dt);
    int n1 = (1100. / _dt);

    init_progressbar(&bar, max_step);
    for (int n=0; n<max_step; n++){

        // square current
        
        if ((n > n0) && (n < n1)){
            for (int i=0; i<info->num_cells * info->cell_ratio[0]; i++){
                cells.ic[i] = 25;
            }
        }

        update_no_delay(n, &cells, &syns, &bck_syns);
        
        write_cell_state(&fdat, &cells);
        write_cell_current(&fdat, &cells);
        write_cell_fire(&fdat, &cells);

        memset(cells.ic, 0, cells.num_cells* sizeof(double));
        progressbar(&bar, n);
    }

    printf("\nDone\n");

    free_cells(&cells);
    free_syns(&syns);
    free_bcksyns(&bck_syns);
}


void init_simulation(simul_info_t *info, neuron_t *cells, syn_t *syns, bcksyn_t *bck_syns)
{
    ntk_t ntk_syns, ntk_bck_syns;

    info->cell_types = (int*) malloc(info->num_cells * sizeof(int));
    gen_cell_type(info->num_cells, info->cell_ratio, info->cell_types);
    

    // init cell structure
    init_cell_vars(cells, info->num_cells, info->cell_types);

    // init syn structures
    init_bi_ntk(info->num_cells, info->num_cells, &ntk_syns); // synaptic network
    gen_bi_random_ntk_with_type(info->cell_types, info->p_syn, info->g_syn, &ntk_syns);
    init_syn_vars(syns, info->num_cells, info->is_delay_on, info->delay, cells->v, &ntk_syns);
    free_bi_ntk(&ntk_syns);

    int *bck_types = (int*) calloc(info->num_bck, sizeof(int));
    init_bi_ntk(info->num_bck, info->num_cells, &ntk_bck_syns); // background synapse
    gen_bi_random_ntk_with_type(bck_types, &(info->p_bck), &(info->g_bck), &ntk_bck_syns);
    init_bcksyn_vars(bck_syns, info->num_cells, info->num_bck, info->tau_bck, &ntk_bck_syns);
    free_bi_ntk(&ntk_bck_syns);

    // set firing probability of the background synapse
    for (int n=0; n<info->num_cells; n++){
        bck_syns->fspk[n] = info->p_fire_bck;
    }
    
    free(bck_types);
}


void empty_simulvar(simul_info_t *info)
{
    for (int n=0; n<_n_types; n++){
        info->cell_ratio[n] = 0;
        info->g_bck[n] = 0;
        info->p_bck[n] = 0;
        info->amp_sq[n] = 0;

        for (int i=0; i<_n_types; i++){
            info->g_syn[n][i] = 0;
            info->p_syn[n][i] = 0;
        }
    }

    info->t_sq_inp[0] = 0;
    info->t_sq_inp[1] = 0;
}


void gen_cell_type(int num_cells, double cell_ratio[], int *cell_types)
{
    int ctp=0, num_last=num_cells*cell_ratio[0];

    for (int n=0; n<num_cells; n++)
    {
        if (n == num_last) { num_last += cell_ratio[++ctp]; }
        cell_types[n] = ctp;
    }
}


void square_current()
{

}


// test tools
void print_ntk_info(ntk_t *ntk, int n_print_node)
{
    for (int n=0; n<n_print_node; n++){
        printf("node %2d, # edges = %d - ", n, ntk->num_edges[n]);
        for (int i=0; i<ntk->num_edges[n]; i++){
            printf("%d ", ntk->adj_list[n][i]);
        }
        printf("\n");
    }
}


void print_post2syn_id(syn_t *syns, int n_print_node)
{
    for (int n=0; n<n_print_node; n++){
        fprintf(stderr, "post id %2d, # syns %2d - ", n, syns->n_post2syn[n]);
        for (int i=0; i<syns->n_post2syn[n]; i++){
            fprintf(stderr, "%d ", syns->id_post2syn[n][i]);
        }
        fprintf(stderr, "\n");
    }
}


void print_variable(double *x, int n_print_node)
{
    for (int n=0; n<n_print_node; n++){
        fprintf(stderr, "%5.2f, ", x[n]);
    }
    fprintf(stderr, "\n");
}

