#include <stdio.h>
#include <math.h>
 

//Random tour generator for all restarts, basen on SPLASH-2 code
//With each thread accessing column wise in the matrix to attcheive
//coalesced access.
__global__ void createToursColumnWise(unsigned short* tourMatrix, 
                            int cities,
                            int restarts,
                            int time){
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
        rand = glo_id + blockIdx.x + time; //blockIdx.x is tourOffset. Check if this is correct
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


__global__ void twoOptKer2(uint32_t* glo_dist, 
                          unsigned short *glo_tours,
                          int* glo_result, 
                          int cities, 
                          int totIter){
    int change, d;
    unsigned int next, i, j, ip1, jp1, block_size, idx;
    block_size = blockDim.x;
    idx = threadIdx.x;
    ChangeTuple localMinChange;
    extern __shared__ unsigned char totShared[];             //shared memory for both tour, minChange and tempRes

    volatile ChangeTuple* tempRes = (volatile ChangeTuple*)&totShared;       //tempRes holds the best local changes found by each thread
    volatile ChangeTuple* minChange = tempRes + block_size;        //minChange holds the current best change 
    volatile unsigned short* tour =
                 (volatile unsigned short*)(minChange + 1);  //tour for this climber
    if(minChange == NULL){
        printf("pointer error\n");
    }

    //copy global tour to shared memory
    for(int t = idx; t < cities+1; t += block_size){
        tour[t] = glo_tours[blockIdx.x * (cities+1) + t];
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
        The i and j index calculated depending on which thread is calculating it.
        ***/
        float tmp;
        for(int ind = idx; ind < totIter; ind += block_size){
            d = 1-(4*(-2*(totIter-ind)));
            //printf("glo: %d, d: %d\n",ind, d);
            tmp = (((-1-(sqrt((float) d)))/2)*(-1))+0.9999;
            //printf("glo: %d, tmp: %f\n",ind, tmp);
            next = (int) tmp;
            //printf("glo: %d, next: %d\n",ind, next);
            i = (cities-2) - (next-1);
            j = (i+2) + (ind-(totIter-((next*(next-1))/2))); //n>>1
            ip1 = i+1;
            jp1 = j+1;
            //printf("glo: %d, i: %d, j: %d\n",ind, i,j);

            change = glo_dist[tour[i]*cities+tour[j]] + 
                    glo_dist[tour[ip1]*cities+tour[jp1]] -
                    (glo_dist[tour[i]*cities+tour[ip1]] +
                    glo_dist[tour[j]*cities+tour[jp1]]);
            //Each thread shall hold the best local change found
            ChangeTuple check = ChangeTuple(change,(unsigned short)i, (unsigned short) j);
            localMinChange = minInd::apply(localMinChange,check);
        }
        //Write each threads local minimum change (best change found)
        //to the shared array tempRes. 
        if(idx < totIter){
            tempRes[idx].change = localMinChange.change;
            tempRes[idx].i = localMinChange.i;
            tempRes[idx].j = localMinChange.j;
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
        while(num_threads != num_elems){
            if (idx < num_threads){
                tempRes[idx] = minInd::apply(tempRes[idx],tempRes[idx + num_threads]);
            }
            __syncthreads();

            num_elems = num_threads;
            num_threads= (num_elems + 1)/ 2;
        }
        ChangeTuple best = minInd::remVolatile(tempRes[0]);
        //Prepare information for swapping
        int temp, swapCities;
        i = best.i + 1;
        j = best.j;
        swapCities = (((j - best.i) + 1) / 2) + i; //the ceiling of j/2 plus i
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

    //copy best local shared memory black to global memory
    for(int t = idx; t < cities+1; t += block_size){
        glo_tours[blockIdx.x * (cities+1) + t] = tour[t];
    }
    
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
