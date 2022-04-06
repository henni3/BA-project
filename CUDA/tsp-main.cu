#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hostSkel.cu.h"
#include "tsp-kernels.cu.h"
#include "dataCollector.cu.h"

#define MAXCITIES 10000


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
    convert<<< num_blocks, block_size >>> (flag_int, flags_d, totIter);
    
    /*
    int* flag = (int*) malloc(totIter*sizeof(int));
    cudaMemcpy(flag, flag_int, totIter*sizeof(int), cudaMemcpyDeviceToHost);
    printf("flag: [");
    for(int i = 0; i < totIter; i++){
        printf("%d, ", flag[i]);
    }
    printf("]\n \n");
    free(flag);

    
    char* flag_c = (char*) malloc(totIter*sizeof(char));
    cudaMemcpy(flag_c, flags_d, totIter*sizeof(char), cudaMemcpyDeviceToHost);
    printf("flag_c: [");
    for(int i = 0; i < totIter; i++){
        printf("%d, ", flag_c[i]);
    }
    printf("]\n");
    free(flag_c);
    */

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
    if (argc != 3) {
        printf("Usage: %s <block-size> <file-name>\n", argv[0]);
        exit(1);
    }
    // Collect input arguments
    int block_size = atoi(argv[1]);
    char* file_name = argv[2];
    
    initHwd();

    // Collect information from datafile into distMatrix and cities
    uint32_t* distMatrix, *kerDist;
    distMatrix = (uint32_t*) malloc(sizeof(uint32_t) * MAXCITIES * MAXCITIES);
    int cities = fileToDistM(file_name, distMatrix);
    if( cities > MAXCITIES){
        printf("too many cities :( \n");
        exit(1);
    }
    distMatrix = (uint32_t*) realloc(distMatrix,sizeof(uint32_t)* cities * cities);
    cudaMalloc((void**)&kerDist, cities*cities*sizeof(uint32_t));
    cudaMemcpy(kerDist, distMatrix, cities*cities*sizeof(uint32_t), cudaMemcpyHostToDevice);

    
    
    /*printf("cities: %d \n", cities );
    printf("matrix: \n");
    for (int i = 0; i < cities; i++){
        for (int j = 0; j < cities; j++){
            printf("%d ", distMatrix[i *cities + j]);
        }
        printf("\n");
    }
    printf("\n");*/


    //Calculate total number of iterations
    int totIter = ((cities-1) * (cities-2))/2;

    //Memory for i-array and j-array
    int *is_d, *js_d;
    cudaMalloc((void**)&is_d, totIter*sizeof(uint32_t));
    cudaMalloc((void**)&js_d, totIter*sizeof(uint32_t));


    init(block_size, cities, totIter, is_d, js_d);

    int* is_h = (int*) malloc(totIter*sizeof(uint32_t));
    cudaMemcpy(is_h, is_d, totIter*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    int k = 0;
    printf("is: [");
    for(int i = 0; i < totIter; i++){
        printf("%d, ", is_h[i]);
        k++;
    }
    printf("]\n");
    printf("k = %d\n", k);

    int* js_h = (int*) malloc(totIter*sizeof(int));
    cudaMemcpy(js_h, js_d, totIter*sizeof(int), cudaMemcpyDeviceToHost);
    printf("js: [");
    for(int i = 0; i < totIter; i++){
        printf("%d, ", js_h[i]);
    }
    printf("]\n");
    free(js_h); free(is_h);
    

    //run 2 opt kernel
    unsigned int num_blocks = (totIter + block_size-1)/block_size; 
    unsigned short *tour, *kerTour;
    tour = (unsigned short*) malloc((cities+1)*sizeof(unsigned short));
    for(int i = 0; i < cities; i++){
        tour[i] = i;
    }

    //tour[0] = 1; tour[1] = 3; tour[2] = 4; tour[3] = 0; tour[4] =2; tour[5] = 1;
    cudaMalloc((void**)&kerTour, (cities+1)*sizeof(unsigned short));
    cudaMemcpy(kerTour, tour, (cities+1)*sizeof(unsigned short), cudaMemcpyHostToDevice);
    unsigned short sharedMemSize = (cities+1) * sizeof(unsigned short) + (block_size*3) * sizeof(int) + 3*sizeof(int);
    twoOptKer<<< num_blocks, block_size, sharedMemSize>>> (kerDist, kerTour, is_d, js_d, cities, totIter);

    

    free(tour); free(distMatrix);
    cudaFree(is_d); cudaFree(js_d);
    cudaFree(kerDist); cudaFree(kerTour);
    return 0;

    
}
