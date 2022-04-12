#include <stdio.h>
      
__global__ void mkIndShp(int* index_shp_d, int len){
    int glb_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(glb_id < len){
        index_shp_d[glb_id] = len - glb_id;
    }
}

__global__ void convert(int* out_flags, char* in_flags, int totIter) {
    int glb_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(glb_id < totIter){
        out_flags[glb_id] = in_flags[glb_id];
    }
}

__global__ void replicate0(int totIter, char* flags_d) {
    int glb_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(glb_id < totIter){
        flags_d[glb_id]=0;
    }
}

__global__ void replicate1(int totIter, int* flags_d) {
    int glb_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(glb_id < totIter){
        flags_d[glb_id]=1;
    }
}

__global__ void mkFlags(int len, int* index_shp_sc_d, char* flags_d) {
    int glb_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(glb_id < len){
        if(glb_id == 0){
             flags_d[0]=1;
        }else{
            flags_d[index_shp_sc_d[glb_id-1]]=1;
        }
    }
}

__global__ void minusOne(int totIter, int* in_arr) {
    int glb_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(glb_id < totIter){
        in_arr[glb_id]=in_arr[glb_id]-1;
    }
}

//Compute the local optimum cost
__device__ int sumTourKernel(uint32_t* glo_dist, 
                                volatile unsigned short *lo_tour, 
                                int cities, 
                                volatile int* result_arr){    
    int idx = threadIdx.x;
    int block_size = blockDim.x;
    int sum = 0;
    int glo_i, glo_ip1; 
    for(int i = idx; i < cities; i += block_size){
        glo_i = lo_tour[i];
        glo_ip1 = lo_tour[i+1];
        sum += glo_dist[glo_i * cities + glo_ip1];
    }
    result_arr[idx] = sum;
    __syncthreads();
    for (int size = block_size /2; size > 0; size /= 2 ){
        if (idx < size) {
            result_arr[idx] += result_arr[idx+size];
        }
        __syncthreads();
    }
    return result_arr[idx];
}

//Random tour generator for all restarts, basen on SPLASH-2 code
__global__ void createTours(unsigned short* tourMatrix, 
                            int cities,
                            int restarts){
    int rand, glo_id, to, temp;
    glo_id = threadIdx.x + blockIdx.x * blockDim.x;
    if(glo_id < restarts){
        //Initiate all tours from 0 to cities
        for(int i = 0; i < cities; i++){
            tourMatrix[(cities+1) * glo_id + i] = i;
        }
        //The last element in all tours is the same as the first (which is 0).
        tourMatrix[(cities+1) * glo_id + cities] = 0;
        
        //Randomize each tour
        rand = glo_id + blockIdx.x; //blockIdx.x is tourOffset. Check if this is correct
        for(int i = 1; i < cities; i++){
            rand = (MULT * rand + ADD) & MASK;
            to = rand % cities;
            if (to <= 0){
                to = 1;
            }
            temp = tourMatrix[(cities+1) * glo_id + i];
            tourMatrix[(cities+1) * glo_id + i] = tourMatrix[(cities+1) * glo_id + to];
            tourMatrix[(cities+1) * glo_id + to] = temp;
        }
        for(int i = 0; i < cities+1; i++){
            printf("glo_id: %d, elem: %d\n",glo_id, tourMatrix[(cities+1) * glo_id + i]);
        }
    }
}


__global__ void twoOptKer(uint32_t* glo_dist, 
                          unsigned short *glo_tour, 
                          int* glo_is, int* glo_js,
                          int* glo_result, 
                          int cities, 
                          int totIter){
    int block_size = blockDim.x;
    int idx = threadIdx.x;
    int glo_id = idx + blockIdx.x * block_size;
    int i, j;
    int32_t localMinChange[3];
    extern __shared__ unsigned char totShared[];             //shared memory for both tour, minChange and tempRes

    volatile int* tempRes = (volatile int*)&totShared;       //tempRes holds the best local changes found by each thread
    volatile int* minChange = tempRes + 3*block_size;        //minChange holds the current best change
    volatile unsigned short* tour =
                 (volatile unsigned short*)(minChange + 3);  //tour for this climber
    if(minChange == NULL){
        printf("pointer error\n");
    }

    /*//Test of shared memory
    int resSize = blockDim.x + cities+1;
    int totSize = resSize+3;
    for(i = threadIdx.x; i < totSize; i += blockDim.x){
        if(i < cities+1){
            tour[i] = glo_tour[i];
            printf("shareTour: %d\n", tour[i]);
        }else if(i > cities && i < resSize){
            int tmp = (i-(cities+1))*3;
            tempRes[tmp] = tmp;
            tempRes[tmp+1] = tmp;  
            tempRes[tmp+2] = tmp;  
            //printf("temp res: fst %d, sec %d, thr %d \n", tempRes[tmp], tempRes[tmp+1], tempRes[tmp+2]);
        }else if(i < totSize){
            int tmpM = (i-resSize)*3;
            minChange[tmpM] = tmpM;
            minChange[tmpM+1] = tmpM;
            minChange[tmpM+2] = tmpM;
            printf("minChange: fst %d, sec %d, thr %d\n", minChange[tmpM], minChange[tmpM+1], minChange[tmpM+2]);
        }
    }*/

    //Preparing data for the 2 opt algorithm
    int ip1, jp1, change;
    //initialize tour to shared memory
    for(int t = idx; t < cities+1; t += block_size){
        tour[t] = glo_tour[t];
    }
    if(idx == 0){
        //initialize minChange to shared memory
        minChange[0] = -1; 
        minChange[1] = 0; 
        minChange[2] = 0;
    }
    
    __syncthreads();
    //Computation for one climber
    while(minChange[0] < 0){
        if(idx < 3){
            minChange[idx] = 0;
        }
        // reset each threads local min change
        localMinChange[0] = 0; 
        localMinChange[1] = 0; 
        localMinChange[2] = 0;
        
        /***
        The 2 opt move
        Each thread calculates the local changes of the given i and j indexes.
        The i and j index are collected (with a stride of block size) from the 
        global i array and in the global j array to acheive coalesecing.
        ***/
        for(int ind = idx; ind < totIter; ind += block_size){
            i = glo_is[ind];
            j = glo_js[ind] + i + 2;
            ip1 = i+1;
            jp1 = j+1;
            change = glo_dist[tour[i]*cities+tour[j]] + 
                    glo_dist[tour[ip1]*cities+tour[jp1]] -
                    (glo_dist[tour[i]*cities+tour[ip1]] +
                    glo_dist[tour[j]*cities+tour[jp1]]);
            //Each thread shall hold the best local change found
            if(change < localMinChange[0]){  
                localMinChange[0] = change; 
                localMinChange[1] = i; 
                localMinChange[2] = j;
            }
        }
        //Write each threads local minimum change (best change found)
        //to the shared array tempRes. 
        if(idx < totIter){
            tempRes[idx*3] = localMinChange[0];
            tempRes[idx*3+1] = localMinChange[1];
            tempRes[idx*3+2] = localMinChange[2];
        }
        __syncthreads();
        
        //Preparation for the reduction on all local minimum changes.
        int num_elems, num_threads;
        if(totIter < block_size){
            num_elems = totIter;
        }else{
            num_elems = block_size;
        }
        num_threads = (num_elems + 1 ) / 2;

        //Reduction on all the local minimum changes found by each thread
        //to find the best minimum change for this climber.
        while(1){
            if(idx < num_threads){
                if (idx + num_threads < num_elems){
                    if (tempRes[idx*3] > tempRes[(idx + num_threads)*3]) {
                        tempRes[idx*3] = tempRes[(idx + num_threads)*3];
                        tempRes[idx*3 + 1] = tempRes[(idx + num_threads)*3 + 1];
                        tempRes[idx*3 + 2] = tempRes[(idx + num_threads)*3 + 2];
                    }
                    else if (tempRes[idx*3] == tempRes[(idx + num_threads)*3]){
                        if (tempRes[idx*3 + 1] > tempRes[(idx + num_threads)*3 +1 ]){
                            tempRes[idx*3] = tempRes[(idx + num_threads)*3 ];
                            tempRes[idx*3 + 1] = tempRes[(idx + num_threads)*3 + 1];
                            tempRes[idx*3 + 2] = tempRes[(idx + num_threads)*3 + 2];
                        }
                        else if (tempRes[idx*3 +1] == tempRes[(idx + num_threads)*3 + 1]){
                            if (tempRes[idx*3 + 2] > tempRes[(idx + num_threads)*3] +2){
                                tempRes[idx*3] = tempRes[(idx + num_threads)*3 ];
                                tempRes[idx*3 + 1] = tempRes[(idx + num_threads)*3 + 1];
                                tempRes[idx*3 + 2] = tempRes[(idx + num_threads)*3 + 2];
                            }
                        }
                    }
                }
            }
            __syncthreads();

            num_elems = num_threads;
            num_threads= (num_elems + 1)/ 2;
            if(num_threads == num_elems){
                break;
            }
        }
        //Prepare information for swapping
        int temp, swapCities;
        i = tempRes[1] + 1;
        j = tempRes[2];
        swapCities = (((tempRes[2] - tempRes[1]) + 1) / 2) + i; //the ceiling of j/2 plus i
        //swap
        for(int t = idx + i; t < swapCities; t += block_size){
            temp = tour[t];
            tour[t] = tour[j - (t - i)];
            tour[j - (t - i)] = temp;
        }
        if(idx < 3){
            minChange[idx] = tempRes[idx];
        }
        __syncthreads();
    }
    
    int local_opt_cost = sumTourKernel(glo_dist, tour, cities, tempRes);
    
    //Writing results to global memory
    if(idx == 0){
        glo_result[blockIdx.x * 2] = local_opt_cost;
        glo_result[blockIdx.x * 2+1] = blockIdx.x;
    }

    
    
    /* Writing results to global memory
    ** THIS DOES NOT WORK AS YOu CAN NOT SYNCRONIZE ACROS GRID.
    if(idx == 0){
        tmp_result[blockIdx.x*2] = local_opt_cost;
        tmp_result[blockIdx.x*2+1] = blockIdx.x;
        atomicAdd(blocksWritten, 1);
    }
    while(blocks_written[0] < num_block);
    if(idx == 0){
        printf("block written: %d, num block: %d", blocks_written, num_block);
    }
    int tot_threads = (num_block + 1 ) / 2;
    int tot_elems = num_block;
    while(1){
        if(glo_id < num_threads){
            if (glo_id + num_threads < num_elems){
                if (tmp_result[glo_id*2] > tmp_result[(glo_id + num_threads)*2]) {
                    tmp_result[glo_id*2] = tmp_result[(glo_id + num_threads)*2];
                    tmp_result[glo_id*2+1] = tmp_result[(glo_id + num_threads)*2+1]
                }
            }
        }
        __syncthreads();

        num_elems = num_threads;
        num_threads= (num_elems + 1)/ 2;
        if(num_threads == num_elems){
            break;
        }
    }
    __syncthreads();
    if(blockIdx.x == tmp_result[1]){
        if(idx == 0){
            glo_result[0] = tmp_result[0];
            glo_result[1] = tour;
        }
    }*/
}

