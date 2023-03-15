#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "hostSkel.cu.h"
#include "tsp-kernels.cu.h"
#include "dataCollector.cu.h"

void init(int block_size, 
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
    mkFlags<<< num_blocks_shp, block_size >>>(len, index_shp_sc_d, flags_d); // was totIter

    cudaMalloc((void**)&flag_int,       totIter*sizeof(int));
    // Convert from char to int.
    convert<<< num_blocks, block_size >>> (flag_int, flags_d, totIter);
    
    //Make i's array
    // 1. scan the flag array
    scanInc<Add<int> > (block_size, totIter, is_d, flag_int, d_tmp_int);
    // 2. minus each element of is_d array with one to get the final is_d array
    minusOne<<< num_blocks, block_size >>> (totIter, is_d);
 
    //Make j's array
    // 1. create an array of ones
    replicate1<<< num_blocks, block_size >>> (totIter, oneArr);
    // 2. segmented scan on the flag array
    sgmScanInc<Add<int> > (block_size, totIter, js_d, flags_d, oneArr, seg_sc_tmp_int, d_tmp_flag);
    // 3. minus each element of js_d array with one to get the final js_d array
    minusOne<<< num_blocks, block_size >>> (totIter, js_d);
    // Store both the information of i's and j's in one array.
    zip<<< num_blocks, block_size>>> (is_d,js_d,totIter);
    
    //free cuda memory
    cudaFree(index_shp_d);  cudaFree(index_shp_sc_d); cudaFree(flag_int);
    cudaFree(flags_d);  cudaFree(d_tmp_int);  cudaFree(d_tmp_flag);
}

void multBlockRed(int *glo_results, int num_blocks_restarts, int block_size, int restarts){
    int num_blocks = (num_blocks_restarts+1)/2;
    int mult_sharedMem = (block_size*2) * sizeof(int);
    int num_elems = restarts;
    while(num_elems > 1){
        multBlockReduce<<<num_blocks, block_size, mult_sharedMem>>>(glo_results, num_elems);
        num_elems = num_blocks;
        num_blocks = (num_blocks + block_size-1)/block_size;
    }
}

void run_original(unsigned short *tourMatrixIn_d, 
                 unsigned short *tourMatrixTrans_d,
                 int *is_d, uint32_t* kerDist, int *glo_results,
                 int block_size, int cities, int restarts, int totIter){
    int num_blocks_restarts, time;
    struct timeval randomTime;

    //Prepare for column wise tour
    num_blocks_restarts = (restarts + block_size-1)/block_size; 
    gettimeofday(&randomTime, NULL);
    time = randomTime.tv_usec;
    //Create tour matrix column wise
    createToursColumnWise<<<num_blocks_restarts, block_size>>> (tourMatrixIn_d, cities, restarts, time);
    transposeTiled<unsigned short, TILE>(tourMatrixIn_d, tourMatrixTrans_d, (cities+1), restarts);
    //compute shared memory size
    size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + 
                            block_size * sizeof(ChangeTuple) + 
                            sizeof(ChangeTuple);
    //run 2 opt kernel 
    twoOptKer<<<restarts, block_size, sharedMemSize>>> (kerDist, tourMatrixTrans_d, 
                                                    is_d, glo_results, 
                                                    cities, totIter);
    //run reduction of all local optimum cost across multiple blocks
    int* glo_res_h = (int*) malloc(block_size*sizeof(int));
    cudaMemcpy(glo_res_h, glo_results, block_size*sizeof(int), cudaMemcpyDeviceToHost);
    printf("glo_res before multBlockRed:  \n[");
    for(int i = 0; i < block_size; i++){
        printf("%d, ", glo_res_h[i]);
    }
    printf("]\n\n");

    multBlockRed(glo_results, num_blocks_restarts, block_size, restarts);
    
    cudaMemcpy(glo_res_h, glo_results, block_size*sizeof(int), cudaMemcpyDeviceToHost);
    printf("glo_res after multBlockRed:  \n[");
    for(int i = 0; i < block_size; i++){
        printf("%d, ", glo_res_h[i]);
    }
    printf("]\n");
    free(glo_res_h); 
}

void run_100cities(unsigned short *tourMatrixIn_d, 
                 unsigned short *tourMatrixTrans_d,
                 int *is_d, uint32_t* kerDist, int *glo_results,
                 int block_size, int cities, int restarts, int totIter){
    int num_blocks_restarts, time;
    struct timeval randomTime;

    //Prepare for column wise tour
    num_blocks_restarts = (restarts + block_size-1)/block_size;
    gettimeofday(&randomTime, NULL);
    time = randomTime.tv_usec;
    //Create tour matrix column wise
    createToursColumnWise<<<num_blocks_restarts, block_size>>> (tourMatrixIn_d, cities, restarts, time);
    transposeTiled<unsigned short, TILE>(tourMatrixIn_d, tourMatrixTrans_d, (cities+1), restarts);
    //Compute shared memory size
    size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + 
                            block_size * sizeof(ChangeTuple) + 
                            sizeof(ChangeTuple) + 
                            cities * cities * sizeof(uint32_t);
    //Run 2 opt kernel
    twoOptKer100Cities<<<restarts, block_size, sharedMemSize>>> (kerDist, tourMatrixTrans_d, 
                                                    is_d, glo_results, 
                                                    cities, totIter);

    //Run reduction of all local optimum cost across multiple blocks
    multBlockRed(glo_results, num_blocks_restarts, block_size, restarts);
}

void run_calculatedIandJ(unsigned short *tourMatrixIn_d, 
                 unsigned short *tourMatrixTrans_d,
                 uint32_t* kerDist, int *glo_results,
                 int block_size, int cities, int restarts, int totIter){
    int num_blocks_restarts, time;
    struct timeval randomTime;

    //Prepare for column wise tour
    num_blocks_restarts = (restarts + block_size-1)/block_size; 
    gettimeofday(&randomTime, NULL);
    time = randomTime.tv_usec;
    //Create tour matrix column wise
    createToursColumnWise<<<num_blocks_restarts, block_size>>> (tourMatrixIn_d, cities, restarts, time);
    transposeTiled<unsigned short, TILE>(tourMatrixIn_d, tourMatrixTrans_d, (cities+1), restarts);
    //Compute shared memory size
    size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + 
                            block_size * sizeof(ChangeTuple) + 
                            sizeof(ChangeTuple);
    //run 2 opt kernel 
    twoOptKerCalculated<<<restarts, block_size, sharedMemSize>>> (kerDist, tourMatrixTrans_d, 
                                                    glo_results, 
                                                    cities, totIter);
    //Run reduction of all local optimum cost across multiple blocks
    multBlockRed(glo_results, num_blocks_restarts, block_size, restarts);
}
