#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hostSkel.cu.h"
#include "tsp-kernels.cu.h"

int init(int block_size, 
         int cities, 
         int totIter, 
         int* is_d, 
         int* js_d){
    int len, *index_shp_d, *index_shp_sc_d, *d_tmp_int;
    int *flag_int, *oneArr, *seg_sc_tmp_int;
    char *flags_d, *d_tmp_flag; 
    
    //Calculate the length of shape array
    len = cities - 2;
    //Calculate block size
    unsigned int num_blocks     = (totIter + block_size-1)/block_size; 
    unsigned int num_blocks_shp = (len + block_size-1)/block_size; 

    //Cuda malloc 
    cudaMalloc((void**)&index_shp_d,    len*sizeof(int));
    cudaMalloc((void**)&index_shp_sc_d, len*sizeof(int));
    cudaMalloc((void**)&flags_d,        totIter*sizeof(char));
    cudaMalloc((void**)&oneArr,         totIter*sizeof(int));

    cudaMalloc((void**)&d_tmp_int,      MAX_BLOCK*sizeof(int));
    cudaMalloc((void**)&seg_sc_tmp_int,   MAX_BLOCK*sizeof(int));
    cudaMalloc((void**)&d_tmp_flag,     totIter*sizeof(char));
    //Create shape array for index
    mkIndShp<<< num_blocks, block_size >>> (index_shp_d, len);

    // Make flag array
    // 1. scan the shape array
    scanInc<Add<int> > (block_size, len, index_shp_sc_d, index_shp_d, d_tmp_int);

    // 2. create an array of zeros
    replicate0<<< num_blocks, block_size >>> (totIter, flags_d);
    
    // 3. scatter the flag array
    mkFlags<<< num_blocks_shp, block_size >>> (totIter, index_shp_sc_d, flags_d);
    cudaMalloc((void**)&flag_int,       totIter*sizeof(int));
    convert<<< num_blocks_shp, block_size >>> (flag_int, flags_d, totIter);
    
    /*
    int* flag = (int*) malloc(totIter*sizeof(int));
    cudaMemcpy(flag, flag_int, totIter*sizeof(int), cudaMemcpyDeviceToHost);
    printf("flag: [");
    for(int i = 0; i < totIter; i++){
        printf("%d, ", flag[i]);
    }
    printf("]\n");
    free(flag); */

    /*
    char* flag_c = (char*) malloc(totIter*sizeof(char));
    cudaMemcpy(flag_c, flags_d, totIter*sizeof(char), cudaMemcpyDeviceToHost);
    printf("flag_c: [");
    for(int i = 0; i < totIter; i++){
        printf("%d, ", flag_c[i]);
    }
    printf("]\n");
    free(flag_c); */


    //Make is array
    // 1. scan the flag array
    scanInc<Add<int> > (block_size, totIter, is_d, flag_int, d_tmp_int);
    // 2. minus each element of is_d array with one to get the final is_d array
     minusOne<<< num_blocks, block_size >>> (totIter, is_d);
 
    //Make js array
    // 1. create an array of ones
    replicate1<<< num_blocks, block_size >>> (totIter, oneArr);
    // 2. segmented scan on the flag array
    sgmScanInc<Add<int> > (block_size, totIter, js_d, flags_d, oneArr, seg_sc_tmp_int, d_tmp_flag);
    // 3. minus each element of js_d array with one to get the final js_d array
    minusOne<<< num_blocks, block_size >>> (totIter, js_d);

    
    //free cuda memory
    cudaFree(index_shp_d);  cudaFree(index_shp_sc_d);
    cudaFree(flags_d);  cudaFree(d_tmp_int);  cudaFree(d_tmp_flag);
    return 0;
}


int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Usage: %s <block-size>\n", argv[0]);
        exit(1);
    }
    
    initHwd();
    int cities = 5;
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
    int *is_d, *js_d;
    cudaMalloc((void**)&is_d, totIter*sizeof(uint32_t));
    cudaMalloc((void**)&js_d, totIter*sizeof(uint32_t));


    init(block_size, cities, totIter, is_d, js_d);
    int* is_h = (int*) malloc(totIter*sizeof(uint32_t));
    cudaMemcpy(is_h, is_d, totIter*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    printf("is: [");
    for(int i = 0; i < totIter; i++){
        printf("%d, ", is_h[i]);
    }
    printf("]\n");

    int* js_h = (int*) malloc(totIter*sizeof(int));
    cudaMemcpy(js_h, js_d, totIter*sizeof(int), cudaMemcpyDeviceToHost);
    printf("js: [");
    for(int i = 0; i < totIter; i++){
        printf("%d, ", js_h[i]);
    }
    printf("]\n");
    free(js_h); free(is_h);

    cudaFree(is_d); cudaFree(js_d);
    return 0;

    
}
