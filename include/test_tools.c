#include <stdio.h>
#include "Izh.h"


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




