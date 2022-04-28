#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "hostSkel.cu.h"
#include "tsp-kernels2.cu.h"
#include "dataCollector.cu.h"

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


    //Prepare for column wise tour
    unsigned short *tourMatrixIn_d, *tourMatrixTrans_d;
    struct timeval randomTime;
    cudaMalloc((void**)&tourMatrixIn_d, (cities+1)*restarts*sizeof(unsigned short));
    cudaMalloc((void**)&tourMatrixTrans_d, (cities+1)*restarts*sizeof(unsigned short));

    unsigned int num_blocks_tour = (restarts + block_size-1)/block_size; 
    gettimeofday(&randomTime, NULL);
    int time = randomTime.tv_usec;

    //Create tour matrix column wise
    createToursColumnWise<<<num_blocks_tour, block_size>>> (tourMatrixIn_d, cities, restarts, time);
    transposeTiled<unsigned short, TILE>(tourMatrixIn_d, tourMatrixTrans_d, (cities+1), restarts);
    cudaFree(tourMatrixIn_d);

    //run 2 opt kernel 
    size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + block_size * sizeof(ChangeTuple) + sizeof(ChangeTuple);
    int *glo_results;
    cudaMalloc((void**)&glo_results, 2*restarts*sizeof(int));
    /*twoOptKer2<<<restarts, block_size, sharedMemSize>>> (kerDist, tourMatrixTrans_d, 
                                                        glo_results, 
                                                        cities, totIter);*/
    //gpuErrchk( cudaPeekAtLastError() );


    //testing timer for twoOptKer2
    int REPEAT;
    int elapsed;
    struct timeval ker2_start, ker2_end, ker2_diff;
    //Dry run
    twoOptKer2<<<restarts, block_size, sharedMemSize>>> (kerDist, tourMatrixTrans_d, 
                                                        glo_results, 
                                                        cities, totIter);
    REPEAT = 9;
    gettimeofday(&ker2_start, NULL); 
    while(REPEAT < 10){
        twoOptKer2<<<restarts, block_size, sharedMemSize>>> (kerDist, tourMatrixTrans_d, 
                                                        glo_results, 
                                                        cities, totIter);
        REPEAT++;
    }
    cudaDeviceSynchronize();
    gettimeofday(&ker2_end, NULL); 
    timeval_subtract(&ker2_diff, &ker2_end, &ker2_start);
    elapsed = (ker2_diff.tv_sec*1e6+ker2_diff.tv_usec) / REPEAT; 
    printf("ker2: Optimized Program runs on GPU in: %lu milisecs, repeats: %d\n", elapsed/1000, REPEAT);
 
    
    //run reduction of all local optimum cost across multiple blocks
    unsigned int num_blocks_gl_re = (num_blocks_tour+1)/2;
    size_t mult_sharedMem = (block_size*2) * sizeof(int);
    for(int i = num_blocks_gl_re; i > 1; i>>=1){
        multBlockReduce<<<i, block_size, mult_sharedMem>>>(glo_results, restarts);
        i++;
    }
    //run reduction on the last block
    multBlockReduce<<<1, block_size, mult_sharedMem>>>(glo_results, restarts);
    cudaDeviceSynchronize();

    //print results
    int* glo_res = (int*) malloc(2*restarts*sizeof(int));
    cudaMemcpy(glo_res, glo_results, 2*restarts*sizeof(int), cudaMemcpyDeviceToHost);
    
    //tour matrix row wise
    unsigned short* tourMatrix_h = (unsigned short*) malloc((cities+1)*restarts*sizeof(unsigned short));
    cudaMemcpy(tourMatrix_h, tourMatrixTrans_d, (cities+1)*restarts*sizeof(unsigned short), cudaMemcpyDeviceToHost);

    
    /*//test tour matrix column wise
    unsigned short* tourMatrixC_h = (unsigned short*) malloc((cities+1)*restarts*sizeof(unsigned short));
    cudaMemcpy(tourMatrixC_h, tourMatrixC_d, (cities+1)*restarts*sizeof(unsigned short), cudaMemcpyDeviceToHost);
    
    printf("Tour R:  [");
    for(int i = 0; i < restarts; i++){
        printf("[");
        for(int j = 0; j < cities+1; j++){
                printf("%d, ", tourMatrix_h[i*(cities+1)+j]);
        }
        printf("]\n");
    }
    printf("]\n\n");

    printf("Tour C:  [");
    for(int i = 0; i < restarts; i++){
        printf("[");
        for(int j = 0; j < cities+1; j++){
                printf("%d, ", tourMatrixC_h[j*restarts+i]);
        }
        printf("]\n");
    }
    printf("]\n");
    free(tourMatrixC_h); cudaFree(tourMatrixC_d);*/

    
    int tourId = glo_res[1];

    printf("Shortest path: %d\n", glo_res[0]);
    printf("Tour:  [");
    for(int i = 0; i < cities+1; i++){
        printf("%d, ", tourMatrix_h[(cities+1)*tourId+i]);
    }
    printf("]\n");


    free(distMatrix); free(tourMatrix_h); free(glo_res); 
    cudaFree(tourMatrixTrans_d);
    cudaFree(kerDist);
    cudaFree(glo_results); 
    printf("done\n");
    return 0;  
}
