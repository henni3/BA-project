#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "hostSkel.cu.h"
#include "tsp-kernels-100.cu.h"
#include "dataCollector.cu.h"

#define CITIES 100

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
    mkFlags<<< num_blocks_shp, block_size >>>(len, index_shp_sc_d, flags_d); // was totIter

    cudaMalloc((void**)&flag_int,       totIter*sizeof(int));
    convert<<< num_blocks, block_size >>> (flag_int, flags_d, totIter);
    
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
    zip<<< num_blocks, block_size>>> (is_d,js_d,totIter);
    
    //free cuda memory
    cudaFree(index_shp_d);  cudaFree(index_shp_sc_d);
    cudaFree(flags_d);  cudaFree(d_tmp_int);  cudaFree(d_tmp_flag);
    return 0;
}

void run_kernels(unsigned short *tourMatrixIn_d, 
                 unsigned short *tourMatrixTrans_d,
                 int *is_d, uint32_t* kerDist, int *glo_results,
                 int block_size, int cities, int restarts, int totIter){
    int num_blocks_tour, num_blocks_gl_re, time;
    struct timeval randomTime;

    //Prepare for column wise tour
    num_blocks_tour = (restarts + block_size-1)/block_size;
    gettimeofday(&randomTime, NULL);
    time = randomTime.tv_usec;
    //Create tour matrix column wise
    createToursColumnWise<<<num_blocks_tour, block_size>>> (tourMatrixIn_d, cities, restarts, time);
    transposeTiled<unsigned short, TILE>(tourMatrixIn_d, tourMatrixTrans_d, (cities+1), restarts);

    //run 2 opt kernel 
    size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + 
                            block_size * sizeof(ChangeTuple) + 
                            sizeof(ChangeTuple) + 
                            cities * cities * sizeof(uint32_t);

    twoOptKer3<<<restarts, block_size, sharedMemSize>>> (kerDist, tourMatrixTrans_d, 
                                                    is_d, glo_results, 
                                                    cities, totIter);

    //run reduction of all local optimum cost across multiple blocks
    num_blocks_gl_re = (num_blocks_tour+1)/2;
    size_t mult_sharedMem = (block_size*2) * sizeof(int);
    for(int i = num_blocks_gl_re; i > 1; i>>=1){
        multBlockReduce<<<i, block_size, mult_sharedMem>>>(glo_results, restarts);
        i++;
    }
    //run reduction on the last block
    multBlockReduce<<<1, block_size, mult_sharedMem>>>(glo_results, restarts);
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Usage: %s <block-size> <file-name> <number-of-restarts>\n", argv[0]);
        exit(1);
    }
    // Collect input arguments
    int block_size = atoi(argv[1]);
    char* file_name = argv[2];
    int restarts = atoi(argv[3]);
    if(restarts <= 0){
        printf("Number of restarts has to be a number larger than 0");
        exit(1);
    }
    
    initHwd();

    //Create varibales
    struct timeval start, end, diff;
    uint32_t* distMatrix, *kerDist;
    int cities, totIter, *is_d, *js_d, *glo_results, *glo_res_h, tourId, elapsed;
    unsigned short *tourMatrixIn_d, *tourMatrixTrans_d, *tourMatrix_h;


    // Collect information from datafile into distMatrix and cities    
    distMatrix = (uint32_t*) malloc(sizeof(uint32_t) * MAXCITIES * MAXCITIES);
    cities = fileToDistM(file_name, distMatrix);
    if( cities > CITIES){
        printf("too many cities :( \n");
        exit(1);
    }
    distMatrix = (uint32_t*) realloc(distMatrix,sizeof(uint32_t) * cities * cities);
    cudaMalloc((void**)&kerDist, cities*cities*sizeof(uint32_t));
    cudaMemcpy(kerDist, distMatrix, cities*cities*sizeof(uint32_t), cudaMemcpyHostToDevice);

    //Calculate total number of iterations
    totIter = ((cities-1) * (cities-2))/2;

    //Cuda malloc
    cudaMalloc((void**)&tourMatrixIn_d, (cities+1)*restarts*sizeof(unsigned short));
    cudaMalloc((void**)&tourMatrixTrans_d, (cities+1)*restarts*sizeof(unsigned short));
    cudaMalloc((void**)&is_d, totIter*sizeof(uint32_t));
    cudaMalloc((void**)&js_d, totIter*sizeof(uint32_t));
    cudaMalloc((void**)&glo_results, 2*restarts*sizeof(int));

    //CPU malloc
    glo_res_h = (int*) malloc(2*sizeof(int));
    tourMatrix_h = (unsigned short*) malloc((cities+1)*restarts*sizeof(unsigned short));

    //testing time for cities 100 program
    gettimeofday(&start, NULL); 
    for(int i = 0; i < GPU_RUNS; i++){
        //run program
        init(block_size, cities, totIter, is_d, js_d);
        run_kernels(tourMatrixIn_d, tourMatrixTrans_d, 
                    is_d, kerDist, glo_results, 
                    block_size, cities, restarts, totIter);

        //get results
        cudaMemcpy(glo_res_h, glo_results, 2*sizeof(int), cudaMemcpyDeviceToHost);
        tourId = glo_res_h[1];
    }

    cudaDeviceSynchronize();
    gettimeofday(&end, NULL); 
    timeval_subtract(&diff, &end, &start);
    elapsed = (diff.tv_sec*1e6+diff.tv_usec) / GPU_RUNS; 
    printf("kernel 100 tour: Optimized Program runs on GPU in: %lu milisecs, repeats: %d\n", elapsed/1000, GPU_RUNS);
    
    //get results
    cudaMemcpy(tourMatrix_h, tourMatrixTrans_d, (cities+1)*restarts*sizeof(unsigned short), cudaMemcpyDeviceToHost);
    
    //print results
    printf("Shortest path: %d\n", glo_res_h[0]);
    printf("Tour:  [");
    for(int i = 0; i < cities+1; i++){
        printf("%d, ", tourMatrix_h[(cities+1)*tourId+i]);
    }
    printf("]\n");
    
    //Clean up
    free(distMatrix); free(tourMatrix_h); free(glo_res_h);  
    cudaFree(is_d); cudaFree(js_d); cudaFree(tourMatrixTrans_d); cudaFree(tourMatrixIn_d);
    cudaFree(kerDist);
    cudaFree(glo_results);
    return 0;
}
