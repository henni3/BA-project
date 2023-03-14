#include "tsp-main-helper.cu.h"

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("Usage: %s <block-size> <file-name> <number-of-restarts> <Which program version? 1 (original), 2 (100Cities) or 3 (calculatedIandJ)>\n", argv[0]);
        exit(1);
    }
    printf("\nBlocksize has to be a number of 2^x otherwise reduce does not work..\n");
    // Collect input arguments
    int block_size = atoi(argv[1]);
    char* file_name = argv[2];
    int restarts = atoi(argv[3]);
    if(restarts <= 0){
        printf("Number of restarts has to be a number larger than 0");
        exit(1);
    }
    int version = atoi(argv[4]);
    if((version < 1) || (version > 3)){
        printf("Wrong program version. You can choose between 1, 2 or 3");
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
    if(((version == 2) && (cities > CITIES)) || (cities > MAXCITIES)){
        printf("too many cities :( \n");
        exit(1);
    }
    distMatrix = (uint32_t*) realloc(distMatrix,sizeof(uint32_t) * cities * cities);
    cudaMalloc((void**)&kerDist, cities*cities*sizeof(uint32_t));
    cudaMemcpy(kerDist, distMatrix, cities*cities*sizeof(uint32_t), cudaMemcpyHostToDevice);
        
    
    //Calculate total number of iterations
    totIter = ((cities-1) * (cities-2))>>1;

    //Cuda malloc
    cudaMalloc((void**)&tourMatrixIn_d, (cities+1)*restarts*sizeof(unsigned short));
    cudaMalloc((void**)&tourMatrixTrans_d, (cities+1)*restarts*sizeof(unsigned short));
    cudaMalloc((void**)&is_d, totIter*sizeof(uint32_t));
    cudaMalloc((void**)&js_d, totIter*sizeof(uint32_t));
    cudaMalloc((void**)&glo_results, 2*restarts*sizeof(int));

    //CPU malloc
    glo_res_h = (int*) malloc(2*sizeof(int));
    tourMatrix_h = (unsigned short*) malloc((cities+1)*restarts*sizeof(unsigned short));
    
    //Dry run init
    init(block_size, cities, totIter, is_d, js_d);
    
    //Testing the original program version
    if(1 == version){
        //Dry run program
        run_original(tourMatrixIn_d, tourMatrixTrans_d, 
                is_d, kerDist, glo_results, 
                block_size, cities, restarts, totIter);
        cudaDeviceSynchronize();
        
        //testing time for original program over GPU runs
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
    
    //Testing the program optimised for 100 cities version
    }else if(2 == version){
        //Dry run program
        run_100cities(tourMatrixIn_d, tourMatrixTrans_d, 
                is_d, kerDist, glo_results, 
                block_size, cities, restarts, totIter);
        cudaDeviceSynchronize();
        
        //testing time for cities 100 program over GPU runs
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

    //Testing the program version where the i and j indexes are calculated (version 3)
    }else{
        //Dry run program
        run_calculatedIandJ(tourMatrixIn_d, tourMatrixTrans_d, 
                kerDist, glo_results, 
                block_size, cities, restarts, totIter);
        cudaDeviceSynchronize();
        
        //testing time for program with calculations of i and j, over GPU runs
        gettimeofday(&start, NULL); 
        for(int i = 0; i < GPU_RUNS; i++){
            //run program
            init(block_size, cities, totIter, is_d, js_d);
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
    
    //Clean up
    free(distMatrix); free(tourMatrix_h); free(glo_res_h);  
    cudaFree(is_d); cudaFree(js_d); cudaFree(tourMatrixTrans_d); cudaFree(tourMatrixIn_d);
    cudaFree(kerDist);
    cudaFree(glo_results);
    return 0;
}
