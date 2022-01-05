#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ntk.h"
#include "mt64.h"


void init_bi_ntk(int num_pre, int num_post, ntk_t *ntk)
{
    ntk->num_pre  = num_pre;
    ntk->num_post = num_post;

    ntk->node_types = (int*) malloc(num_pre * sizeof(int));
    ntk->num_edges  = (int*) calloc(num_pre, sizeof(int));

    ntk->adj_list = (int**) malloc(num_pre * sizeof(int*));
    ntk->strength = (double**) malloc(num_pre * sizeof(double*));

    for (int n=0; n<num_pre; n++){
        ntk->adj_list[n] = (int*) malloc(_ntk_block_size * sizeof(int));
        ntk->strength[n] = (double*) malloc(_ntk_block_size * sizeof(double));
    }
}


void gen_bi_random_ntk_with_type(int *pre_node_types, int *post_node_types,
                            double p_cnt_type[_n_types][_n_types],
                            double str_cnt_type[_n_types][_n_types], ntk_t *ntk)
{
    // node_type, 1 ~ ntype
    // p_cnt_type, ntype x ntype
    int i, j, pre_tp, post_tp; // i: pre id; j: post id
    double p, strength;

    memcpy(ntk->node_types, pre_node_types, ntk->num_pre * sizeof(int));

    for (i=0; i<ntk->num_pre; i++){
        pre_tp = pre_node_types[i];
        for (j=0; j<ntk->num_post; j++){
            post_tp = post_node_types[j];

            // printf("%d:%d - %d:%d\n", ntk->num_pre, pre_tp, ntk->num_post, post_tp);

            if ((ntk->num_pre == ntk->num_post) && (i == j)){
                continue;
            }

            p = genrand64_real2();
            if (p < p_cnt_type[pre_tp][post_tp]){
                // connection
                append_node(j, ntk->num_edges+i, ntk->adj_list+i);
                ntk->num_edges[i]--;
                
                // append coupling strength
                strength = str_cnt_type[pre_tp][post_tp];
                append_value(strength, ntk->num_edges+i, ntk->strength+i);
            }
        }
    }
}


void gen_bi_random_ntk_fixed_indeg(int *pre_node_types, int *post_node_types,
                                double p_cnt_type[_n_types][_n_types],
                                double str_cnt_type[_n_types][_n_types], ntk_t *ntk)
{
    memcpy(ntk->node_types, pre_node_types, ntk->num_pre * sizeof(int));
    
    int num_types[_n_types] = {0,};
    int id_boundary[_n_types] = {0,};
    int *ptr_ctp = pre_node_types;

    for (int i=0; i<ntk->num_pre; i++){
        id_boundary[*ptr_ctp] = i;
        num_types[*ptr_ctp]++; 
        ptr_ctp++;        
    }

    for (int i=_n_types-1; i>0; i--){
        id_boundary[i] = id_boundary[i-1]+1;
    }
    id_boundary[0] = 0;

    int num_in_deg;
    int post_tp;
    int *used = (int*) malloc(ntk->num_pre * sizeof(int));

    for (int i=0; i<ntk->num_post; i++){
        post_tp = post_node_types[i];
        for (int pre_tp=0; pre_tp<_n_types; pre_tp++){
            num_in_deg = p_cnt_type[pre_tp][post_tp] * num_types[pre_tp];
            memset(used, 0, ntk->num_pre * sizeof(int));

            int pre_id, j=0, a=id_boundary[pre_tp], num=num_types[pre_tp];
            double s = str_cnt_type[pre_tp][post_tp];

            while (j < num_in_deg){
                pre_id = genrand64_real2()*num + a;
                if (used[pre_id] == 0){

                    int *ptr_num_edge = ntk->num_edges + pre_id;

                    append_node(i, ptr_num_edge, ntk->adj_list+pre_id);
                    (*ptr_num_edge)--;
                    append_value(s, ptr_num_edge, ntk->strength+pre_id);

                    used[pre_id] = 1;
                    j++;
                }
            }
        }
    }
}


void free_bi_ntk(ntk_t *ntk)
{   
    free(ntk->node_types);
    free(ntk->num_edges);
    
    for (int n=0; n<ntk->num_pre; n++){
        free(ntk->adj_list[n]);
        free(ntk->strength[n]);
    }
    free(ntk->adj_list);
    free(ntk->strength);
}


void append_node(int node_id, int *n_edge, int **adj_list)
{
    // adj_list: pointer of the n th adj list
    if ((*n_edge > 0) && (*n_edge % _ntk_block_size == 0)){ // increase the size
        *adj_list = (int*) realloc(*adj_list, *n_edge + _ntk_block_size);
    }
    (*adj_list)[*n_edge] = node_id;
    (*n_edge)++;
}


void append_value(double val, int *n_edge, double **val_list)
{
    // adj_list: pointer of the n th adj list
    if ((*n_edge > 0) && ((*n_edge)%_ntk_block_size == 0)){ // increase the size
        *val_list = (double*) realloc(*val_list, *n_edge + _ntk_block_size);
    }
    (*val_list)[*n_edge] = val;
    (*n_edge)++;
}