#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <errno.h>
#include "../lib/lib_Izh.h"
#include "../lib/mt64.h"
#include "../lib/parson.h"

extern int errno;
typedef double fourDouble[4];

// FILE IO obs
FILE *fid_v, *fid_u, *fid_r, *fid_id;

int main(int argc, char **argv){

    // init with random seed
    init_genrand64(100);

    // set simulation args
    double dt=0.01, tmax=100;
    int N,i,j,n,nitr;

    char *finfo, *fout;

    // read args
    int c;
    while ((c=getopt(argc, argv, "f:o:t:")) != -1){
        switch (c) {
            case 'f':
                finfo = optarg;
                break;
            case 'o':
                fout = optarg;
                break;
            case 't':
                tmax = atof(optarg);
                break;
        }
    }

    // printf("f input = %s, f output = %s\n", finfo, fout);

    // read args

}

// json obj
void read_simulation_args(char fjson[], int *N, int *Ntypes, fourDouble **CellParams){
    printf("Info fname = %s\n", fjson);

    int i,j;
    JSON_Value *rootValue;
    JSON_Object *rootObject;
    JSON_Array *json_arr;

    if ((rootValue = json_parse_file(fjson)) == NULL){
        errno = ENOENT;
        perror("error");
        exit(0);
    }

    rootObject = json_value_get_object(rootValue);

    // read the # of cells
    *N = json_array_get_number(rootObject, "NumCells");
    *Ntypes = json_array_get_number(rootObject, "NumTypes");

    *CellParams = (fourDouble*) malloc(sizeof(fourDouble)*(*Ntypes));
    for (i=0; i<*Ntypes; i++){
        json_arr = json_object_get_array(rootObject, "CellParams");
        for (j=0; j<4; j++){
            (*CellParams)[i][j] = json_array_get_number(json_arr, j);
        }
    }
    

}
