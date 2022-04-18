#ifndef _NTK
#define _NTK

#define _ntk_block_size 500
#define MAX_TYPE 6

typedef struct _ntk_t {

    int num_pre, num_post;
    int *node_types;
    int *num_edges;
    int **adj_list;
    double **strength;

} ntk_t;


void init_bi_ntk(int num_pre, int num_post, ntk_t *ntk);
void gen_bi_random_ntk_mean_deg(int *pre_node_types, int *post_node_types, 
                            int mean_out_deg[MAX_TYPE][MAX_TYPE],
                            double str_cnt_type[MAX_TYPE][MAX_TYPE], ntk_t *ntk);
void gen_bi_random_ntk_fixed_indeg(int *pre_node_types, int *post_node_types,
                            int mean_out_deg[MAX_TYPE][MAX_TYPE],
                            double str_cnt_type[MAX_TYPE][MAX_TYPE], ntk_t *ntk);
void gen_bi_random_ntk_with_type(int *pre_node_types, int *post_node_types,
                            double p_cnt_type[MAX_TYPE][MAX_TYPE],
                            double str_cnt_type[MAX_TYPE][MAX_TYPE], ntk_t *ntk);
void free_bi_ntk(ntk_t *ntk);
void append_node(int node_id, int *n_edge, int **adj_list);
void append_value(double val, int *n_edge, double **val_list);


#endif
