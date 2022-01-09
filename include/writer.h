#ifndef WRITER
#define WRITER


typedef enum _WRITER_VAR
{
    V_ONLY = 1,
    U_ONLY = 2,
    I_ONLY = 4,
    SPK_ONLY = 8

} WRITER_VAR;


typedef struct _writer_t
{
    FILE *v;
    FILE *u;
    FILE *i;
    FILE *t_spk;

} writer_t;




#endif