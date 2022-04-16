#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hostSkel.cu.h"
#include "tsp-kernels.cu.h"
#include "dataCollector.cu.h"

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
    //cudaDeviceSynchronize();
    /*int* indSha = (int*) malloc(len*sizeof(int));
    cudaMemcpy(indSha, index_shp_d, len*sizeof(int), cudaMemcpyDeviceToHost);
    printf("indSha: [");
    for(int i = 0; i < len; i++){
        printf("%d, ", indSha[i]);
    }
    printf("]\n \n");
    free(indSha);*/
    // Make flag array
    // 1. scan the shape array
    scanInc<Add<int> > (block_size, len, index_shp_sc_d, index_shp_d, d_tmp_int);
    //gpuErrchk( cudaPeekAtLastError() );
    //cudaDeviceSynchronize();
    /*int* scan = (int*) malloc(len*sizeof(int));
    cudaMemcpy(scan, index_shp_sc_d, len*sizeof(int), cudaMemcpyDeviceToHost);
    printf("scan: [");
    for(int i = 0; i < len; i++){
        printf("%d, ", scan[i]);
    }
    printf("]\n \n");
    free(scan);*/

    // 2. create an array of zeros
    replicate0<<< num_blocks, block_size >>> (totIter, flags_d);
    
    // 3. scatter the flag array
    mkFlags<<< num_blocks_shp, block_size >>>(len, index_shp_sc_d, flags_d); // was totIter
    //gpuErrchk( cudaPeekAtLastError() );
    //gpuErrchk( cudaDeviceSynchronize() );
    cudaMalloc((void**)&flag_int,       totIter*sizeof(int));
    convert<<< num_blocks, block_size >>> (flag_int, flags_d, totIter);
    
    
    /*int* flag = (int*) malloc(totIter*sizeof(int));
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
    free(flag_c);*/

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
    if (argc != 4) {
        printf("Usage: %s <block-size> <file-name> <number-of-restarts>\n", argv[0]);
        exit(1);
    }
    // Collect input arguments
    int block_size = atoi(argv[1]);
    //printf("block size %d, \n", block_size);
    char* file_name = argv[2];
    int restarts = atoi(argv[3]);
    if(restarts <= 0){
        printf("Number of restarts has to be a number larger than 0");
        exit(1);
    }
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

    /*int* is_h = (int*) malloc(totIter*sizeof(uint32_t));
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
    free(js_h); free(is_h);*/
    

    //Create tour matrix
    unsigned short *tourMatrix;
    cudaMalloc((void**)&tourMatrix, (cities+1)*restarts*sizeof(unsigned short));
    unsigned int num_blocks_tour = (restarts + block_size-1)/block_size; 
    createTours<<<num_blocks_tour, block_size>>> (tourMatrix, cities, restarts);

    size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + (block_size*3) * sizeof(int) + 3*sizeof(int);
    printf("sharedmemSize used in twoOptKer : %d \n", sharedMemSize);

    int *glo_results;
    cudaMalloc((void**)&glo_results, 2*restarts*sizeof(int));

    //run 2 opt kernel 
    twoOptKer<<<restarts, block_size, sharedMemSize>>> (kerDist, tourMatrix, is_d, js_d, glo_results, cities, totIter);
    //gpuErrchk( cudaPeekAtLastError() );
    printf("after twoOptKernel\n");
        
    size_t gl_re_sharedMem = (block_size*2) * sizeof(int);
    unsigned int num_blocks_gl_re = (num_blocks_tour+1)/2;
    //run reduction of all local optimum cost
    size_t mult_sharedMem = (block_size*2) * sizeof(int);
    for(num_blocks_gl_re; num_blocks_gl_re > 1; num_blocks_gl_re>>=1){
        multBlockReduce<<<num_blocks_gl_re, block_size, gl_re_sharedMem>>>(glo_results, restarts);
        num_blocks_gl_re++;
    }
    multBlockReduce<<<1, block_size, gl_re_sharedMem>>>(glo_results, restarts);
    
    int* glo_res = (int*) malloc(2*restarts*sizeof(int));
    cudaMemcpy(glo_res, glo_results, 2*restarts*sizeof(int), cudaMemcpyDeviceToHost);
    printf("results:\n");
   

    printf("result: %d, block_id\n", glo_res[0], glo_res[1]);
    free(glo_res);

    free(distMatrix);
    cudaFree(is_d); cudaFree(js_d); cudaFree(tourMatrix);
    cudaFree(kerDist);
    cudaFree(glo_results); 
    return 0;

    
}
