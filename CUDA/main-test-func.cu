#include "tsp-main-helper.cu.h"

int main() {
    initHwd();

    int block_size, cities, totIter, restarts;
    int *is_d, *js_d, *is_h;
    unsigned short *tourMatrixIn_d, *tourMatrixIn_h;
    block_size = 1024;
    cities = 5;
    restarts = 5;
    totIter = ((cities-1) * (cities-2))/2;

    //TEST: is i and j array correct? Yes it is
    cudaMalloc((void**)&is_d, totIter*sizeof(uint32_t));
    cudaMalloc((void**)&js_d, totIter*sizeof(uint32_t));
    is_h = (int*) malloc(totIter*sizeof(uint32_t));

    init(block_size, cities, totIter, is_d, js_d);

    /*cudaMemcpy(is_h, is_d, totIter*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    printf("i     j\n");
    for(int ind = 0; ind < totIter; ind++){
        int num = is_h[ind];
        i = num >> 16;
        j = (num & 0xffff) + i + 2;
        printf("%d     %d\n", i,j);
    }
    printf("end\n");*/

    
    //TEST: Does createTours work correctly? Nope
    cudaMalloc((void**)&tourMatrixIn_d, (cities+1)*restarts*sizeof(unsigned short));
    tourMatrixIn_h = (unsigned short*) malloc((cities+1)*restarts*sizeof(unsigned short));

    int num_blocks_tour = (restarts + block_size-1)/block_size; 
    //struct timeval randomTime;
    //gettimeofday(&randomTime, NULL);
    //int time = randomTime.tv_usec;
    int time = 4000;
    createToursColumnWise<<<num_blocks_tour, block_size>>> (tourMatrixIn_d, cities, restarts, time);
    
    cudaMemcpy(tourMatrixIn_h, tourMatrixIn_d, (cities+1)*restarts*sizeof(unsigned short), cudaMemcpyDeviceToHost);
    printf("Tour matrix\n");
    for(int i = 0; i < cities+1; i++){
        for(int j = 0; j < restarts; j++){
            printf("%hu, ", tourMatrixIn_h[i * restarts + j]);
        }
        printf("\n");
    }
    printf("end");


    cudaFree(is_d); cudaFree(js_d); cudaFree(tourMatrixIn_d);
    free(is_h); free(tourMatrixIn_h);

    return 0;
}