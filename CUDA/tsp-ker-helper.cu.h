#include <math.h>

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

__global__ void zip(int* array1, int* array2, int size){
    int glb_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(glb_id < size){
        array1[glb_id] = (array1[glb_id] << 16) | array2[glb_id];
    }
}


//Compute the local optimum cost
__device__ int sumTourKernel(uint32_t* glo_dist, 
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
    }

    result_arr[idx].change = sum;
    __syncthreads();
    for (int size = block_size >> 1; size > 0; size >>= 1 ){
        if (idx < size) {
            result_arr[idx].change += result_arr[idx+size].change;

        }
        __syncthreads();
    }
    return result_arr[0].change;
}

//Random tour generator for all restarts, basen on SPLASH-2 code
//With each thread accessing column wise in the matrix to attcheive
//coalesced access.
__global__ void createToursColumnWise(unsigned short* tourMatrix,
                            int cities,
                            int restarts,
                            int time){
    uint32_t rand, glo_id, temp, to;
    glo_id = threadIdx.x + blockIdx.x * blockDim.x;
    if(glo_id < restarts){
        //Initiate all tours from 0 to cities
        for(int i = 0; i < cities; i++){
            tourMatrix[restarts * i + glo_id] = i;
        }
        //The last city in all tours is the same as the first (which is city 0).
        tourMatrix[restarts * cities + glo_id] = 0;
        
        //Randomize each tour
        rand = glo_id + blockIdx.x + time;
        for(int i = 1; i < cities; i++){
            rand = (MULT * rand + ADD) & MASK;
            to = rand % cities;
            if(to <= 0){
                to = 1;
            }
            //toMat[restarts * i + glo_id] = to;
            temp = tourMatrix[restarts * i + glo_id];
            tourMatrix[restarts * i + glo_id] = tourMatrix[restarts * to + glo_id];
            tourMatrix[restarts * to + glo_id] = temp;
        }
    }
}


//Reduction on all the local minimum changes found by each thread
//to find the best minimum change for this climber.
__device__ void reduceLocalMinChange(int block_size, 
                                    volatile ChangeTuple* tempRes){
    int idx = threadIdx.x;
    for (int size = block_size >> 1; size > 0; size >>= 1 ){
        if(idx < size){
            tempRes[idx] = minInd::apply(tempRes[idx],tempRes[idx + size]);
        }
        __syncthreads();
    }
}


/** Reduces all the local optimum cost to find the
 *  global optimum cost. 
 *  This is done by parallel reduction across multiple blocks. 
 * **/
__global__ void multBlockReduce(int* glo_result, 
                                int num_elems){
    int idx, glo_id, block_size, elem1, elem2;
    idx = threadIdx.x;
    block_size = blockDim.x;
    glo_id = idx + (block_size * 2) * blockIdx.x;
    extern __shared__ int sharedMem[];    //shared memory
    
    //Padding with max integer value
    /*if(glo_id > num_elems){
        glo_result[glo_id*2] = INT_MAX;
        glo_result[(glo_id + tot_threads)*2] = INT_MAX;
    }else if((glo_id + tot_threads) > num_elems){
        glo_result[(glo_id + tot_threads)*2] = INT_MAX;
    }*/

    /*if(glo_id > num_elems){
        elem1 = INT_MAX;
        elem2 = INT_MAX;
    }else if((glo_id + block_size) > num_elems){
        elem1 = glo_result[glo_id*2];
        elem2 = INT_MAX;
    }else{
        elem1 = glo_result[glo_id*2];
        elem2 = glo_result[(glo_id + block_size)*2];
    }*/
    
    if(glo_id > (num_elems-1)){
        sharedMem[idx*2] = INT_MAX;
        sharedMem[(idx*2)+1] = INT_MAX;
    }else if((glo_id + block_size) > (num_elems-1)){
        sharedMem[idx*2] = glo_result[(glo_id*2)];
        sharedMem[(idx*2)+1] = glo_result[(glo_id*2)+1];
    }else{
        elem1 = glo_result[glo_id*2];
        elem2 = glo_result[(glo_id + block_size)*2];
        if(elem1 <= elem2){
            sharedMem[idx*2] = elem1;
            sharedMem[(idx*2)+1] = glo_result[(glo_id*2)+1];
        }else{
            sharedMem[idx*2] = elem2;
            sharedMem[(idx*2)+1] = glo_result[((glo_id + block_size)*2)+1];
        }
    }
    __syncthreads();

    glo_result[(glo_id*2)] = sharedMem[idx*2];
    glo_result[(glo_id*2)+1] = sharedMem[(idx*2)+1];

    //reduce on elements in shared memory
    /*for (int size = block_size >> 1; size > 0; size >>= 1 ){
        if(idx < size){
            if(sharedMem[idx*2] > sharedMem[(idx + size)*2]){
                sharedMem[idx*2] = sharedMem[(idx + size)*2];
                sharedMem[(idx*2)+1] = sharedMem[((idx + size)*2)+1];
            }
        }
        __syncthreads();
    }
    if(idx == 0){        
        glo_result[blockIdx.x*2] = sharedMem[0];
        glo_result[(blockIdx.x*2)+1] = sharedMem[1];
    }*/

    /*
    //Find limit of how many elements are to be reduced by this block.
    int tot_threads, n;
    if(num_elems < ((block_size* 2) * (blockIdx.x + 1))){
        n = num_elems - (block_size * 2) * blockIdx.x;
    }else{
        n = block_size * 2;
    }
    tot_threads = (n + 1) >> 1;

    extern __shared__ int sharedMem[];    //shared memory
    //Compare element in global memory across two blocks and the smallest
    //element will be written to shared memory (first reduce layer).
    if(idx < tot_threads){
        elem1 = glo_result[glo_id*2];
        if(idx + tot_threads < n){
            elem2 = glo_result[(glo_id + tot_threads)*2];
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
            if(idx + i < n){
                if(sharedMem[idx*2] > sharedMem[(idx + i)*2]){
                    sharedMem[idx*2] = sharedMem[(idx + i)*2];
                    sharedMem[(idx*2)+1] = sharedMem[((idx + i)*2)+1];
                }
            }
            n = i;
            i++;
            __syncthreads();
        }
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
    }*/
}

