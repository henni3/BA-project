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

__device__ int sumTourKernel2(uint32_t* glo_dist, 
                                volatile unsigned short *lo_tour, 
                                int cities, 
                                volatile ChangeTuple* result_arr){    
    int idx = threadIdx.x;
    int block_size = blockDim.x;
    int sum = 0;
    int glo_i, glo_ip1; 
    for(int i = idx; i < cities; i += block_size){
        glo_i = lo_tour[i];
        glo_ip1 = lo_tour[i+1];
        sum += glo_dist[glo_i * cities + glo_ip1];
        if(blockIdx.x < 1){
            printf("sum is %d with i %d and i+1 %d \n", sum, glo_i, glo_ip1 );
        }
    }

    result_arr[idx].change = sum;
    if (blockIdx.x < 1){
        printf("resultarr after init %d, with thread %d \n", result_arr[idx].change, idx);

    }
    __syncthreads();
    for (int size = block_size /2; size > 0; size /= 2 ){
        if (idx < size) {
            result_arr[idx].change += result_arr[idx+size].change;
            //printf("rsult arr %d \n", result_arr[idx].change );
        }
        __syncthreads();
    }
    if (idx < 1) {
        //printf("computed result is %d \n", result_arr[idx].change);
    }
    return result_arr[idx].change;
}

//Random tour generator for all restarts, basen on SPLASH-2 code
//With each thread accessing row wise in the matrix. This does not
//attcheive coalesced access.
__global__ void createToursRowWise(unsigned short* tourMatrix, 
                            int cities,
                            int restarts,
                            int time){
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
        rand = glo_id + blockIdx.x + time; //blockIdx.x is tourOffset. Check if this is correct
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
    }
}


//Random tour generator for all restarts, basen on SPLASH-2 code
//With each thread accessing column wise in the matrix to attcheive
//coalesced access.
//!!!!!!!!!!!!!!!!!!!!!!!!!!NOOOOT DONE!!!!!!!!!!!!!!!!!!
__global__ void createToursColumnWise(unsigned short* tourMatrix, 
                            int cities,
                            int restarts){
    int rand, glo_id, to, temp;
    glo_id = threadIdx.x + blockIdx.x * blockDim.x;
    if(glo_id < restarts){
        //Initiate all tours from 0 to cities
        for(int i = 0; i < cities; i++){
            tourMatrix[restarts * i + glo_id] = i;
        }
        //The last element in all tours is the same as the first (which is 0).
        tourMatrix[restarts * cities + glo_id] = 0;
        
        //Randomize each tour
        rand = glo_id + blockIdx.x; //blockIdx.x is tourOffset. Check if this is correct
        for(int i = 1; i < cities; i++){
            rand = (MULT * rand + ADD) & MASK;
            to = rand % cities;
            if (to <= 0){
                to = 1;
            }
            temp = tourMatrix[restarts * i + glo_id];
            tourMatrix[restarts * i + glo_id] = tourMatrix[restarts * to + glo_id];
            tourMatrix[restarts * to + glo_id] = temp;
        }
        //TRANSPOSE MATRIX!!
    }
}

__global__ void twoOptKer2(uint32_t* glo_dist, 
                          unsigned short *glo_tours, 
                          int* glo_is, int* glo_js,
                          int* glo_result, 
                          int cities, 
                          int totIter){
    int block_size = blockDim.x;
    int idx = threadIdx.x;
    int i, j;
    ChangeTuple localMinChange;
    extern __shared__ unsigned char totShared[];             //shared memory for both tour, minChange and tempRes

    volatile ChangeTuple* tempRes = (volatile ChangeTuple*)&totShared;       //tempRes holds the best local changes found by each thread
    volatile ChangeTuple* minChange = tempRes + block_size;        //minChange holds the current best change (cos?)
    volatile unsigned short* tour =
                 (volatile unsigned short*)(minChange);  //tour for this climber
    if(minChange == NULL){
        printf("pointer error\n");
    }

    //Preparing data for the 2 opt algorithm
    int ip1, jp1, change;
    //copy global tour to shared memory
    for(int t = idx; t < cities+1; t += block_size){
        tour[t] = glo_tours[blockIdx.x * (cities+1) + t];
        printf("tour is %d, with idx %d \n", tour[t], idx);
    }
    if(idx == 0){
        //initialize minChange to shared memory
        minChange[0] = ChangeTuple();
        minChange[0].change = -1;
    }
    
    __syncthreads();
    //Computation for one climber
    while(minChange[0].change < 0){
        if(idx < 1){
            minChange[0] = ChangeTuple();
        }
        // reset each threads local min change
        localMinChange = ChangeTuple();
        
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
            if (blockIdx.x < 1) {
                    printf(" change is %d \n", change);
                }
            //Each thread shall hold the best local change found
            if(change < localMinChange.change){
                if (blockIdx.x < 1) {
                    printf("best change found %d \n", change);
                }
                localMinChange.change = change;
                localMinChange.i = i;
                localMinChange.j = j;
            }
        }
        //Write each threads local minimum change (best change found)
        //to the shared array tempRes. 
        if(idx < totIter){
            tempRes[idx].change = localMinChange.change;
            tempRes[idx].i = localMinChange.i;
            tempRes[idx].j = localMinChange.j;
            if (blockIdx.x < 1) {
                    printf("tempre change :D %d, with idx %d\n", tempRes[idx].change, idx);
                }
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
            /*if(idx < num_threads){
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
            }*/
            //printf("we get here \n");
            if (idx < num_threads){
                tempRes[idx] = minInd::apply(tempRes[idx],tempRes[idx + num_threads]);
                if (blockIdx.x < 1) {
                    printf("tempres change %d tempres i %d tempres j %d \n", tempRes[idx].change,tempRes[idx].i, tempRes[idx].j );
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
        i = tempRes[0].i + 1;
        j = tempRes[0].j;
        swapCities = (((j - tempRes[0].i) + 1) / 2) + i; //the ceiling of j/2 plus i
        //swap
        for(int t = idx + i; t < swapCities; t += block_size){
            temp = tour[t];
            tour[t] = tour[j - (t - i)];
            tour[j - (t - i)] = temp;
        }
        if(idx < 1){
            minChange[idx].change = tempRes[idx].change;
            minChange[idx].j = tempRes[idx].j;
            minChange[idx].i = tempRes[idx].i;
        }
        __syncthreads();
    }
    
    int local_opt_cost = sumTourKernel2(glo_dist, tour, cities, tempRes);
    if(idx < 1){
        printf("local opt cost %d \n", local_opt_cost);
    }
    
    
    //Writing local optimum results to global memory
    if(idx == 0){
        glo_result[blockIdx.x * 2] = local_opt_cost;
        glo_result[blockIdx.x * 2+1] = blockIdx.x;
    }
}



__global__ void twoOptKer(uint32_t* glo_dist, 
                          unsigned short *glo_tours, 
                          int* glo_is, int* glo_js,
                          int* glo_result, 
                          int cities, 
                          int totIter){
    int block_size = blockDim.x;
    int idx = threadIdx.x;
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

    //Preparing data for the 2 opt algorithm
    int ip1, jp1, change;
    //copy global tour to shared memory
    for(int t = idx; t < cities+1; t += block_size){
        tour[t] = glo_tours[blockIdx.x * (cities+1) + t];
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
    
    //Writing local optimum results to global memory
    if(idx == 0){
        glo_result[blockIdx.x * 2] = local_opt_cost;
        glo_result[blockIdx.x * 2+1] = blockIdx.x;
    }
}





/** Reduces all the local optimum cost to find the
 *  global optimum cost. 
 *  This is done by parallel reduction across multiple blocks. **/

__global__ void multBlockReduce(int* glo_result, 
                                int num_elem){
    int n, idx, glo_id, tot_threads;
    idx = threadIdx.x;
    glo_id = idx + (blockDim.x * 2) * blockIdx.x;
    
    //Find limit of how many elements are to be reduced by this block.
    if(num_elem < ((blockDim.x * 2) * (blockIdx.x + 1))){
        n = num_elem - (blockDim.x * 2) * blockIdx.x;
    }else{
        n = blockDim.x * 2;
    }
    tot_threads = (n + 1) / 2;

    extern __shared__ int sharedMem[];    //shared memory
    //Compare element in global memory across two blocks and the smallest
    //element will be written to shared memory (first reduce layer).
    if(idx < tot_threads){
        int elem1 = glo_result[glo_id*2];
        if(idx + tot_threads < n){
            int elem2 = glo_result[(glo_id + tot_threads)*2];
            if(elem1 <= elem2){
                sharedMem[idx*2] = elem1;
                sharedMem[(idx*2)+1] = glo_result[(glo_id*2)+1];
            }else{
                sharedMem[idx*2] = elem2;
                sharedMem[(idx*2)+1] = glo_result[((glo_id + tot_threads)*2)+1];
            }
        }else{
            sharedMem[idx*2] = elem1;
            sharedMem[(idx*2)+1] = glo_result[(glo_id*2)+1];
        }
        __syncthreads();
        n = tot_threads;
        tot_threads = (n+1)/2;

        //reduce on elements in shared memory (second layer to second 
        //last layer of reduce).
        for(int i = tot_threads; i > 1; i>>=1){
            if(idx < i){
                if(idx + i < n){
                    if(sharedMem[idx*2] > sharedMem[(idx + i)*2]){
                        sharedMem[idx*2] = sharedMem[(idx + i)*2];
                        sharedMem[(idx*2)+1] = sharedMem[((idx + i)*2)+1];
                    }
                }
            }
            n = i;
            i++;
            __syncthreads();
        }
        //__syncthreads();
        //Compare the last two elements of the last reduce layer and
        //write to global memory.
        if(idx == 0){
            if(sharedMem[0] > sharedMem[2]){
                glo_result[blockIdx.x*2] = sharedMem[2];
                glo_result[(blockIdx.x*2)+1] = sharedMem[3];
            }else{
                glo_result[blockIdx.x*2] = sharedMem[0];
                glo_result[(blockIdx.x*2)+1] = sharedMem[1];
            }
        }
    }
}
