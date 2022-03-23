#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hostSkel.cu.h"
#include "tsp-kernels.cu.h"

int init(int block_size, int cities, int totIter, int* is, int* js){
    int len, *index_shp_d, *index_shp_sc_d, *d_tmp_int;
    char  *flags_d, *d_tmp_flag; 

    //Calculate the length of shape array
    len = cities - 2;
    //initHwd();
    //Calculate block size
    unsigned int num_blocks     = (totIter + block_size-1)/block_size; 
    unsigned int num_blocks_shp = (len + block_size-1)/block_size; 

    //Cuda malloc 
    cudaMalloc((void**)&index_shp_d,    len*sizeof(int));
    cudaMalloc((void**)&index_shp_sc_d, len*sizeof(int));
    cudaMalloc((void**)&flags_d,        totIter*sizeof(char));

    cudaMalloc((void**)&d_tmp_int,   MAX_BLOCK*sizeof(int));
    cudaMalloc((void**)&d_tmp_flag,  MAX_BLOCK*sizeof(char));
    
    //Create shape array for index
    mkIndShp<<< num_blocks, block_size >>> (index_shp_d, len);
    printf("hej-main11\n");
    // Make flag array
    // 1. scan the shape array
    scanInc<Add<int> > (block_size, len, index_shp_sc_d, index_shp_d, d_tmp_int);
    printf("hej\n");
    // 2. create an array of zeros
    replicate0<<< num_blocks, block_size >>> (totIter, flags_d);

    // 3. scatter the flag array
    mkFlags<<< num_blocks_shp, block_size >>> (len, index_shp_sc_d, flags_d);
    printf("index_shape_scan: [");
    for (int i = 0; i < len; i++) {
        printf("%d, ", index_shp_sc_d[i]);
    }
    printf("] \n"); 

    printf("index_shape: [");
    for (int i = 0; i < len; i++) {
        printf("%d, ", index_shp_d[i]);
    }
    printf("] \n");

    printf("flag_Arr: [");
    for (int i = 0; i < len; i++) {
        printf("%d, ", flags_d[i]);
    }
    printf("] \n");

    //Make is array
    // 1. scan the shape array
    scanInc<Add<int> > (block_size, len, index_shp_sc_d, index_shp_d, d_tmp_int);


    //free cuda memory
    cudaFree(index_shp_d);  cudaFree(index_shp_sc_d);
    cudaFree(flags_d);  cudaFree(d_tmp_int);  cudaFree(d_tmp_flag);
    return 0;
}


int main(int argc, char* argv[]) {
    initHwd();
    printf("hej");
    int cities = 330;
    //Calculate total number of iterations
    int totIter = ((cities-1) * (cities-2))/2;

    int block_size = atoi(argv[1]);
   // uint32_t totDist = cities * cities;
    //uint32_t* distM = (uint32_t*) malloc((totDist)*sizeof(uint32_t));
    //uint32_t* tour = (uint32_t*) malloc((cities + 1 ) * sizeof(uint32_t));
    //uint32_t tempDist[25] = {0,4,6,8,3,4,0,4,5,2,6,4,0,2,3,8,5,2,0,4,3,2,3,4,0};
    //uint32_t tempTour[6] = {4,2,0,3,1,4};

    //memcpy(distM, tempDist, sizeof(uint32_t) * (totDist));
    //memcpy(tour, tempTour, sizeof(uint32_t) * (cities+1));
    
    //Memory for i-array and j-array
    int *is, *js;
    cudaMalloc((void**)&is, totIter*sizeof(uint32_t));
    cudaMalloc((void**)&js, totIter*sizeof(uint32_t));
    
    init(block_size, cities, totIter, is, js);

    free(is); free(js);
    return 0;

    
}
