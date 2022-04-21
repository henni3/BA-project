


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
    }

    result_arr[idx].change = sum;
    __syncthreads();
    for (int size = block_size /2; size > 0; size /= 2 ){
        if (idx < size) {
            result_arr[idx].change += result_arr[idx+size].change;

        }
        __syncthreads();
    }
    return result_arr[idx].change;
}



