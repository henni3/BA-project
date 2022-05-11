#include "tsp-main-helper.cu.h"

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
    struct timeval randomTime, start, end, diff;
    uint32_t* distMatrix, *kerDist, num_blocks_tour, num_blocks_gl_re;
    int cities, totIter, *is_d, *js_d, *glo_results, *glo_res_h, tourId, REPEAT, elapsed;
    unsigned short *tourMatrixIn_d, *tourMatrixTrans_d, *tourMatrix_h;
    size_t mult_sharedMem;


    // Collect information from datafile into distMatrix and cities    
    distMatrix = (uint32_t*) malloc(sizeof(uint32_t) * MAXCITIES * MAXCITIES);
    cities = fileToDistM(file_name, distMatrix);
    if( cities > MAXCITIES){
        printf("too many cities :( \n");
        exit(1);
    }
    distMatrix = (uint32_t*) realloc(distMatrix,sizeof(uint32_t)* cities * cities);
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
    int* restart_array;
    cudaMalloc((void**)&restart_array, restarts * sizeof(int));

    //CPU malloc
    glo_res_h = (int*) malloc(2*restarts*sizeof(int));
    tourMatrix_h = (unsigned short*) malloc((cities+1)*restarts*sizeof(unsigned short));

    //testing timer for cities 100 program
    REPEAT = 0;
    gettimeofday(&start, NULL); 
    while(REPEAT < 10){
        init(block_size, cities, totIter, is_d, js_d);

        //Prepare for column wise tour
        num_blocks_tour = (restarts + block_size-1)/block_size; 
        gettimeofday(&randomTime, NULL);
        int time = randomTime.tv_usec * 10;
        //Create tour matrix column wise
        createToursColumnWise<<<num_blocks_tour, block_size>>> (tourMatrixIn_d, cities, restarts, time);
        transposeTiled<unsigned short, TILE>(tourMatrixIn_d, tourMatrixTrans_d, (cities+1), restarts);
        //printf("size of change tuple = %d \n", sizeof(ChangeTuple));
        //run 2 opt kernel 
        size_t sharedMemSize = (cities+1) * sizeof(unsigned short) + block_size * sizeof(ChangeTuple) + sizeof(ChangeTuple);
        //printf("sharedmemSize used in twoOptKer : %d \n", sharedMemSize);

        twoOptKer<<<restarts, block_size, sharedMemSize>>> (kerDist, tourMatrixTrans_d, 
                                                        is_d, glo_results, 
                                                        cities, totIter, restart_array);
        //run reduction of all local optimum cost across multiple blocks
        num_blocks_gl_re = (num_blocks_tour+1)/2;
        mult_sharedMem = (block_size*2) * sizeof(int);
        for(int i = num_blocks_gl_re; i > 1; i>>=1){
            multBlockReduce<<<i, block_size, mult_sharedMem>>>(glo_results, restarts);
            i++;
        }
        //run reduction on the last block
        multBlockReduce<<<1, block_size, mult_sharedMem>>>(glo_results, restarts);

        //print results
        cudaMemcpy(glo_res_h, glo_results, 2*restarts*sizeof(int), cudaMemcpyDeviceToHost);
        int* host_restart = (int*) malloc(restarts * sizeof(int));
        cudaMemcpy(host_restart, restart_array, restarts* sizeof(int),cudaMemcpyDeviceToHost);
        int re_sum = 0;
        for (int i = 0; i < restarts; i++){
            re_sum += host_restart[i];
        }
        float average = (float) re_sum /  (float) restarts;
        printf("average nr. of restarts is %f, for %d climbers \n", average, restarts);
        cudaFree(restart_array);
        free(host_restart);

        
        //tour matrix row wise
        cudaMemcpy(tourMatrix_h, tourMatrixTrans_d, (cities+1)*restarts*sizeof(unsigned short), cudaMemcpyDeviceToHost);
        
        tourId = glo_res_h[1];
        REPEAT++;
    }
    cudaDeviceSynchronize();
    gettimeofday(&end, NULL); 
    timeval_subtract(&diff, &end, &start);
    elapsed = (diff.tv_sec*1e6+diff.tv_usec) / REPEAT; 
    printf("Original kernel: Optimized Program runs on GPU in: %lu milisecs, repeats: %d\n", elapsed/1000, REPEAT);
    
    printf("Shortest path: %d\n", glo_res_h[0]);
    printf("Tour:  [");
    for(int i = 0; i < cities+1; i++){
        printf("%d, ", tourMatrix_h[(cities+1)*tourId+i]);
    }
    printf("]\n");
    
    cudaFree(tourMatrixIn_d);
    free(distMatrix); free(tourMatrix_h); free(glo_res_h); 
    cudaFree(is_d); cudaFree(js_d); cudaFree(tourMatrixTrans_d); 
    cudaFree(kerDist);
    cudaFree(glo_results);
    return 0;
}