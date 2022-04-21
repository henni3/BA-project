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
    

    //Create tour matrix row wise
    unsigned short *tourMatrixR_d; //*tourMatrixC_d;
    struct timeval randomTime;
    cudaMalloc((void**)&tourMatrixR_d, (cities+1)*restarts*sizeof(unsigned short));
    unsigned int num_blocks_tour = (restarts + block_size-1)/block_size; 
    gettimeofday(&randomTime, NULL);
    int time = randomTime.tv_usec;
    createToursRowWise<<<num_blocks_tour, block_size>>> (tourMatrixR_d, cities, restarts, time);

    /*//Create tour matrix column wise
    cudaMalloc((void**)&tourMatrixC_d, (cities+1)*restarts*sizeof(unsigned short));
    createToursColumnWise<<<num_blocks_tour, block_size>>> (tourMatrixC_d, cities, restarts);*/

    //run 2 opt kernel 
    //size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + (block_size*3) * sizeof(int) + 3*sizeof(int);
    size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + block_size * sizeof(ChangeTuple) + sizeof(ChangeTuple);
    //printf("sharedmemSize used in twoOptKer : %d \n", sharedMemSize);
    int *glo_results;
    cudaMalloc((void**)&glo_results, 2*restarts*sizeof(int));
    twoOptKer2<<<restarts, block_size, sharedMemSize>>> (kerDist, tourMatrixR_d, 
                                                        is_d, js_d, glo_results, 
                                                        cities, totIter);
    //gpuErrchk( cudaPeekAtLastError() );
 
    
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
    cudaMemcpy(tourMatrix_h, tourMatrixR_d, (cities+1)*restarts*sizeof(unsigned short), cudaMemcpyDeviceToHost);

    
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
    cudaFree(tourMatrixR_d);
    cudaFree(kerDist);
    cudaFree(glo_results); 
    return 0;

    
}
