#include "tsp-main-helper.cu.h"

int main() {
    initHwd();

    //TEST: is i and j array correct?
    int block_size, cities, totIter, *is_d, *js_d, *is_h;
    block_size = 32;
    cities = 5;
    totIter = ((cities-1) * (cities-2))/2;

    cudaMalloc((void**)&is_d, totIter*sizeof(uint32_t));
    cudaMalloc((void**)&js_d, totIter*sizeof(uint32_t));
    is_h = (int*) malloc(totIter*sizeof(uint32_t));

    init(block_size, cities, totIter, is_d, js_d);

    cudaMemcpy(is_h, is_d, totIter*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    printf("i     j\n")
    for(int ind = 0; ind < totIter; ind++){
        int num = is_h[ind];
        i = num >> 16;
        j = (num & 0xffff) + i + 2;
        printf("%d     %d\n", i,j);
    }
    printf("end\n")

    return 0;
}