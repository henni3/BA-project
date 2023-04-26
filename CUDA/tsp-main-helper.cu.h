#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "hostSkel.cu.h"
#include "tsp-kernels.cu.h"
#include "dataCollector.cu.h"

/************************************************************** 
 *  Init() computes the i and j indexes and stores them in 
 *  is_d and js_d.
 **************************************************************/
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

/************************************************************** 
 *  Reduces all the local optimum cost to find the global
 *  optimum cost. 
 *  This is done by parallel reduction across multiple blocks. 
 **************************************************************/
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

/******************************************************* 
 * run_original runs the 2-opt program with no 
 * optimizations.
*******************************************************/
void run_original(unsigned short *tourMatrixIn_d, 
                 unsigned short *tourMatrixTrans_d,
                 int *is_d, uint32_t* kerDist, int *glo_results,
                 int block_size, int cities, int restarts, int totIter){
  
    int num_blocks_restarts; 
    num_blocks_restarts = (restarts + block_size-1)/block_size; // nvm
    //Create randomized tours
    createToursColumnWise<<<num_blocks_restarts, block_size>>> (tourMatrixIn_d, cities, restarts);
    transposeTiled<unsigned short, TILE>(tourMatrixIn_d, tourMatrixTrans_d, (cities+1), restarts);

    //compute shared memory size
    size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + 
                            block_size * sizeof(ChangeTuple) + 
                            sizeof(ChangeTuple);
    //run 2 opt kernel 
    twoOptKer<<<restarts, block_size, sharedMemSize>>> (kerDist, 
                                                    tourMatrixTrans_d, 
                                                    is_d, glo_results, 
                                                    cities, totIter);

    /*int* host_tmp = (int*)malloc(restarts*2*sizeof(int));
    cudaMemcpy(host_tmp, glo_results, restarts*2*sizeof(int), cudaMemcpyDeviceToHost);
    printf("\n\nPath lengths before reduce:\n");    
    for(int i=0; i<restarts; i++) {
        printf("%d, ", host_tmp[2*i]);
    }
    printf("\n\n");
    */
    //run reduction of all local optimum cost across multiple blocks
    multBlockRed(glo_results, num_blocks_restarts, block_size, restarts);
}

/******************************************************* 
 * run_100cities runs the 2-opt program with max 100 cities
 * where tours is stored in shared memory in the kernel.
*******************************************************/
void run_100cities(unsigned short *tourMatrixIn_d, 
                 unsigned short *tourMatrixTrans_d,
                 int *is_d, uint32_t* kerDist, int *glo_results,
                 int block_size, int cities, int restarts, int totIter){
    int num_blocks_restarts;
    num_blocks_restarts = (restarts + block_size-1)/block_size;
    //Create randomized tours
    createToursColumnWise<<<num_blocks_restarts, block_size>>> (tourMatrixIn_d, cities, restarts);
    transposeTiled<unsigned short, TILE>(tourMatrixIn_d, tourMatrixTrans_d, (cities+1), restarts);
    
    //Compute shared memory size
    size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + 
                            block_size * sizeof(ChangeTuple) + 
                            sizeof(ChangeTuple) + 
                            cities * cities * sizeof(uint32_t);
    //Run 2 opt kernel
    twoOptKer100Cities<<<restarts, block_size, sharedMemSize>>> (kerDist, 
                                                    tourMatrixTrans_d, 
                                                    is_d, glo_results, 
                                                    cities, totIter);

    //Run reduction of all local optimum cost across multiple blocks
    multBlockRed(glo_results, num_blocks_restarts, block_size, restarts);
}

/******************************************************* 
 * run_calculatedIandJ() runs the 2-opt program where
 * index i and index j is calculated
*******************************************************/
void run_calculatedIandJ(unsigned short *tourMatrixIn_d, 
                 unsigned short *tourMatrixTrans_d,
                 uint32_t* kerDist, int *glo_results,
                 int block_size, int cities, int restarts, int totIter){
    int num_blocks_restarts;
    num_blocks_restarts = (restarts + block_size-1)/block_size;
    
    //Create randomized tours
    createToursColumnWise<<<num_blocks_restarts, block_size>>> (tourMatrixIn_d, cities, restarts);
    transposeTiled<unsigned short, TILE>(tourMatrixIn_d, tourMatrixTrans_d, (cities+1), restarts);
    
    //Compute shared memory size
    size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + 
                            block_size * sizeof(ChangeTuple) + 
                            sizeof(ChangeTuple);
    //run 2 opt kernel 
    twoOptKerCalculated<<<restarts, block_size, sharedMemSize>>> (kerDist, 
                                                    tourMatrixTrans_d,
                                                    glo_results, 
                                                    cities, totIter);
    //Run reduction of all local optimum cost across multiple blocks
    multBlockRed(glo_results, num_blocks_restarts, block_size, restarts);
}

/******************************************************* 
 * test runs the 2-opt program with no 
 * optimizations.
*******************************************************/
void run_test(unsigned short *tourMatrixIn_d, 
                 unsigned short *tourMatrixTrans_d,
                 int *is_d, uint32_t* kerDist, int *glo_results, int *counter, 
                 int block_size, int cities, int restarts, int totIter){
  
    int num_blocks_restarts; 
    num_blocks_restarts = (restarts + block_size-1)/block_size; // nvm
    //Create randomized tours
    createToursColumnWise<<<num_blocks_restarts, block_size>>> (tourMatrixIn_d, cities, restarts);
    transposeTiled<unsigned short, TILE>(tourMatrixIn_d, tourMatrixTrans_d, (cities+1), restarts);

    //compute shared memory size
    size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + 
                            block_size * sizeof(ChangeTuple) + 
                            sizeof(ChangeTuple);
    //run 2 opt kernel 
    twoOptKer_test<<<restarts, block_size, sharedMemSize>>> (kerDist, 
                                                    tourMatrixTrans_d, 
                                                    is_d, glo_results, counter, 
                                                    cities, totIter);

    /*int* host_tmp = (int*)malloc(restarts*2*sizeof(int));
    cudaMemcpy(host_tmp, glo_results, restarts*2*sizeof(int), cudaMemcpyDeviceToHost);
    printf("\n\nPath lengths before reduce:\n");    
    for(int i=0; i<restarts; i++) {
        printf("%d, ", host_tmp[2*i]);
    }
    printf("\n\n");
    */
    //run reduction of all local optimum cost across multiple blocks
    multBlockRed(glo_results, num_blocks_restarts, block_size, restarts);
    multBlockRed(counter, num_blocks_restarts, block_size, restarts);
}


/******************************************************* 
 * runProgram() prepares variables and calls functions
 * and kernels to compute the shortest rout in TSP using
 * 2-opt algorithm. 
*******************************************************/
void runProgram(char* file_name, int restarts, int version){
    //Create varibales
    struct timeval start, end, diff;
    uint32_t* distMatrix, *kerDist;
    int block_size, cities, totIter, *is_d, *js_d, *glo_results, *glo_res_h, tourId, elapsed;
    unsigned short *tourMatrixIn_d, *tourMatrixTrans_d, *tourMatrix_h;

    // for testing 
    int *counter, *counter_h;

    // Collect information from datafile into distMatrix and cities    
    distMatrix = (uint32_t*) malloc(sizeof(uint32_t) * MAXCITIES * MAXCITIES);
    cities = fileToDistM(file_name, distMatrix);
    if(((version == 2) && (cities > CITIES)) || (cities > MAXCITIES)){
        printf("too many cities :( \n");
        exit(1);
    }
    distMatrix = (uint32_t*) realloc(distMatrix,sizeof(uint32_t) * cities * cities);
    cudaMalloc((void**)&kerDist, cities*cities*sizeof(uint32_t));
    cudaMemcpy(kerDist, distMatrix, cities*cities*sizeof(uint32_t), cudaMemcpyHostToDevice);
        
    
    //Calculate total number of iterations
    totIter = ((cities-1) * (cities-2))>>1;
    //Calculate block size for kernels.
    if(totIter > 512){
        block_size = 1024;
    }else if(totIter > 256){
        block_size = 512;
    }else if(totIter > 128){
        block_size = 256;
    }else if(totIter > 64){
        block_size = 128;
    }else if(totIter > 32){
        block_size = 64;
    }else{
        block_size = 32;
    }

    //Cuda malloc
    cudaMalloc((void**)&tourMatrixIn_d, (cities+1)*restarts*sizeof(unsigned short));
    cudaMalloc((void**)&tourMatrixTrans_d, (cities+1)*restarts*sizeof(unsigned short));
    cudaMalloc((void**)&is_d, totIter*sizeof(uint32_t));
    cudaMalloc((void**)&js_d, totIter*sizeof(uint32_t));
    cudaMalloc((void**)&glo_results, 2*restarts*sizeof(int));
    //testing
    cudaMalloc(&counter, 2 * restarts * sizeof(int));


    //CPU malloc
    glo_res_h = (int*) malloc(2*sizeof(int));
    tourMatrix_h = (unsigned short*) malloc((cities+1)*restarts*sizeof(unsigned short));

    // testing

    counter_h = (int*) malloc(2*sizeof(int));
    
    //Dry run init
    init(block_size, cities, totIter, is_d, js_d);
    
    //Running the original program version
    if(1 == version){
        //Dry run program
        run_original(tourMatrixIn_d, tourMatrixTrans_d, 
                is_d, kerDist, glo_results, 
                block_size, cities, restarts, totIter);
        cudaDeviceSynchronize();
        
        //Taking time for original program over GPU runs
        gettimeofday(&start, NULL); 
        for(int i = 0; i < GPU_RUNS; i++){
            //run program
            init(block_size, cities, totIter, is_d, js_d);
            run_original(tourMatrixIn_d, tourMatrixTrans_d, 
                        is_d, kerDist, glo_results, 
                        block_size, cities, restarts, totIter);

            //get results
            cudaMemcpy(glo_res_h, glo_results, 2*sizeof(int), cudaMemcpyDeviceToHost);
            tourId = glo_res_h[1];
        }
        cudaDeviceSynchronize();
        gettimeofday(&end, NULL);
    
    //Running the program optimised for 100 cities version
    }else if(2 == version){
        //Dry run program
        run_100cities(tourMatrixIn_d, tourMatrixTrans_d, 
                is_d, kerDist, glo_results, 
                block_size, cities, restarts, totIter);
        cudaDeviceSynchronize();
        
        //Taking time for cities 100 program over GPU runs
        gettimeofday(&start, NULL); 
        for(int i = 0; i < GPU_RUNS; i++){
            //run program
            init(block_size, cities, totIter, is_d, js_d);
            run_100cities(tourMatrixIn_d, tourMatrixTrans_d, 
                        is_d, kerDist, glo_results, 
                        block_size, cities, restarts, totIter);

            //get results
            cudaMemcpy(glo_res_h, glo_results, 2*sizeof(int), cudaMemcpyDeviceToHost);
            tourId = glo_res_h[1];
        }
        cudaDeviceSynchronize();
        gettimeofday(&end, NULL);

    //Running the program version where the i and j indexes are calculated (version 3)
    }else if( 3 == version){
        //Dry run program
        run_calculatedIandJ(tourMatrixIn_d, tourMatrixTrans_d, 
                kerDist, glo_results, 
                block_size, cities, restarts, totIter);
        cudaDeviceSynchronize();
        
        //Taking time for program with calculations of i and j, over GPU runs
        gettimeofday(&start, NULL); 
        for(int i = 0; i < GPU_RUNS; i++){
            //run program
            run_calculatedIandJ(tourMatrixIn_d, tourMatrixTrans_d, 
                        kerDist, glo_results, 
                        block_size, cities, restarts, totIter);

            //get results
            cudaMemcpy(glo_res_h, glo_results, 2*sizeof(int), cudaMemcpyDeviceToHost);
            tourId = glo_res_h[1];
        }
        cudaDeviceSynchronize();
        gettimeofday(&end, NULL);
    }
    else if ( 4 == version) {
                //Dry run program
        run_test(tourMatrixIn_d, tourMatrixTrans_d, 
                is_d, kerDist, glo_results, counter,
                block_size, cities, restarts, totIter);
        cudaDeviceSynchronize();
        
        //Taking time for original program over GPU runs
        gettimeofday(&start, NULL); 
        for(int i = 0; i < GPU_RUNS; i++){
            //run program
            init(block_size, cities, totIter, is_d, js_d);
            run_original(tourMatrixIn_d, tourMatrixTrans_d, 
                        is_d, kerDist, glo_results, 
                        block_size, cities, restarts, totIter);

            //get results
            cudaMemcpy(glo_res_h, glo_results, 2*sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(counter_h, counter, 2 * sizeof(int), cudaMemcpyDeviceToHost);
            tourId = glo_res_h[1];
        }
        cudaDeviceSynchronize();
        gettimeofday(&end, NULL);

    }
    else {
        printf("version not recoginized, something went wrong, understood version: %d", version);
        exit(-1);
    }
    
    
 
    timeval_subtract(&diff, &end, &start);
    elapsed = (diff.tv_sec*1e6+diff.tv_usec) / GPU_RUNS; 
    printf("Version: %d. Optimized program runs on GPU in: %lu milisecs, repeats: %d\n", version, elapsed/1000, GPU_RUNS);
    
    //get results
    cudaMemcpy(tourMatrix_h, tourMatrixTrans_d, (cities+1)*restarts*sizeof(unsigned short), cudaMemcpyDeviceToHost);
    
    //print results
    printf("Shortest path: %d\n", glo_res_h[0]);
    printf("Tour:  [");
    for(int i = 0; i < cities+1; i++){
        printf("%d, ", tourMatrix_h[(cities+1)*tourId+i]);
    }
    printf("]\n");

    if ( version == 4) {
        int while_tot = counter_h[0];
        printf("number of while iteartions across all blocks = %d \n", while_tot);
    }
    
    //Clean up
    free(distMatrix); free(tourMatrix_h); free(glo_res_h); free(counter_h);
    cudaFree(is_d); cudaFree(js_d); cudaFree(tourMatrixTrans_d); cudaFree(tourMatrixIn_d);
    cudaFree(kerDist);
    cudaFree(glo_results);
    cudaFree(counter);  
}