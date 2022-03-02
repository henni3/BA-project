#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hostSkel.cu.h"
#include "kernels.cu.h"

int twoOptMove(int block_size, int cities){
    int   totIter, *index_shp_d, *index_shp_sc_d, *d_tmp_int;
    char  *flags_d, *d_tmp_flag; 

    //Calculate total number of iterations
    totIter = ((cities-1) * (cities-2))/2;

    //Calculate the length of shape array
    len = cities - 2;

    //Cuda malloc 
    cudaMalloc((void**)&index_shp_d,    len*sizeof(int));
    cudaMalloc((void**)&index_shp_sc_d, len*sizeof(int));
    cudaMalloc((void**)&flags_d,        totIter*sizeof(char));

    cudaMalloc((void**)&d_tmp_int,   MAX_BLOCK*sizeof(int));
    cudaMalloc((void**)&d_tmp_flag,  MAX_BLOCK*sizeof(char));

    //Create shape array for index
    for(int i = 0; i < len; i++){
        index_shp_d[i] = len - i;
    }

    // Make flag array
    // 1. scan the shape array
    scanInc< Add<int> > (block_size, totIter, index_shp_sc_d, index_shp_d, d_tmp_int);

    // 2. create an array of zeros
    replicate0<<< num_blocks, block_size >>> ( totIter, flags_d );

    // 3. scatter the flag array
    mkFlags<<< num_blocks_shp, block_size >>> ( mat_rows, mat_shp_sc_d, flags_d );

    //free cuda memory
    cudaFree(index_shp_d);  cudaFree(index_shp_sc_d);
    cudaFree(flags_d);  cudaFree(d_tmp_int);  cudaFree(d_tmp_flag);
}


int main(int argc, char* argv[]) {
    int cities = 5;
    uint32_t totDist = cities * cities;
    uint32_t* distM = (uint32_t*) malloc((totDist)*sizeof(uint32_t));
    uint32_t* tour = (uint32_t*) malloc((cities + 1 ) * sizeof(uint32_t));

    memcpy(distM, (uint32_t[25]){0,4,6,8,3,4,0,4,5,2,6,4,0,2,3,8,5,2,0,4,3,2,3,4,0}, sizeof(uint32_t) * (totDist));
    memcpy(tour, (int[6]) {4,2,0,3,1,4}, sizeof(int) * (cities+1));

    
}