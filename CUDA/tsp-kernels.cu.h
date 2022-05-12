#include "tsp-ker-helper.cu.h"

__global__ void twoOptKer(uint32_t* glo_dist, 
                          unsigned short *glo_tours, 
                          int* glo_is,
                          int* glo_result, 
                          int cities, 
                          int totIter){
    int block_size = blockDim.x;
    int idx = threadIdx.x;
    int i, j, change;
    int repeats = 0;
    ChangeTuple localMinChange;
    extern __shared__ unsigned char totShared[];             //shared memory for both tour, minChange and tempRes

    volatile ChangeTuple* tempRes = (volatile ChangeTuple*)&totShared;       //tempRes holds the best local changes found by each thread
    volatile ChangeTuple* minChange = tempRes + block_size;        //minChange holds the current best change
    volatile unsigned short* tour =
                 (volatile unsigned short*)(minChange + 1);  //tour for this climber
    if(minChange == NULL){
        printf("pointer error\n");
    }

    //Preparing data for the 2 opt algorithm
    int ip1, jp1;

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
           //repeats++;
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
        //float tmp;
        for(int ind = idx; ind < totIter; ind += block_size){
            int num = glo_is[ind];
            i = num >> 16;
            j = (num & 0xffff) + i + 2;
            ip1 = i+1;
            jp1 = j+1; 
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
        //scanIncBlock<minInd>((typename minInd::RedElTP*) tempRes,(unsigned int) idx);
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
    
    int local_opt_cost = sumTourKernel(glo_dist, tour, cities, tempRes);


    //copy best local shared memory black to global memory
    for(int t = idx; t < cities+1; t += block_size){
        glo_tours[blockIdx.x * (cities+1) + t] = tour[t];
    }
    
    //Writing local optimum results to global memory
    if(idx == 0){
        glo_result[blockIdx.x * 2] = local_opt_cost;
        glo_result[blockIdx.x * 2+1] = blockIdx.x;
        //re_array[blockIdx.x] = repeats;
    }
}




__global__ void twoOptKer100Cities(uint32_t* glo_dist, 
                          unsigned short *glo_tours, 
                          int* glo_is,
                          int* glo_result, 
                          int cities, 
                          int totIter){
    int block_size = blockDim.x;
    int idx = threadIdx.x;
    int i, j, change, ip1, jp1;

    ChangeTuple localMinChange;
    extern __shared__ unsigned char totShared[];             //shared memory for both tour, minChange and tempRes
    volatile ChangeTuple* tempRes = (volatile ChangeTuple*)&totShared;       //tempRes holds the best local changes found by each thread
    volatile ChangeTuple* minChange = tempRes + block_size; //minChange holds the current best change (cos?)
    volatile uint32_t* shared_Dist = (volatile uint32_t*) (minChange + 1);       
    volatile unsigned short* tour =
                 (volatile unsigned short*)(shared_Dist + cities * cities);  //tour for this climber
    if(minChange == NULL){
        printf("pointer error\n");
    }

    ChangeTuple tmpMinChange;

    //copy gloabal dist to shared memory
    for (int t = idx; t < cities * cities; t += block_size) {
        shared_Dist[t] = glo_dist[t];
    }

    //copy global tour to shared memory
    for(int t = idx; t < cities+1; t += block_size){
        tour[t] = glo_tours[blockIdx.x * (cities+1) + t];
    }

    //initialize minChange to shared memory
    tmpMinChange.change = -1;

    __syncthreads();

    //Computation for one climber
    while(tmpMinChange.change < 0){
        tmpMinChange = ChangeTuple();

        // reset each threads local min change
        localMinChange = ChangeTuple();
        
        /***
        The 2 opt move
        Each thread calculates the local changes of the given i and j indexes.
        The i and j index are collected (with a stride of block size) from the 
        global i array and in the global j array to acheive coalesecing.
        ***/
        for(int ind = idx; ind < totIter; ind += block_size){
            int num = glo_is[ind];
            i = num >> 16;
            j = (num & 0xffff) + i + 2;
            ip1 = i+1;
            jp1 = j+1;
            change = shared_Dist[tour[i]*cities+tour[j]] + 
                    shared_Dist[tour[ip1]*cities+tour[jp1]] -
                    (shared_Dist[tour[i]*cities+tour[ip1]] +
                    shared_Dist[tour[j]*cities+tour[jp1]]);

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
        
        tmpMinChange.change = tempRes[idx].change;
        tmpMinChange.j = tempRes[idx].j;
        tmpMinChange.i = tempRes[idx].i;
        __syncthreads();
    }
    
    int local_opt_cost = sumTourKernel(glo_dist, tour, cities, tempRes);

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



__global__ void twoOptKerCalculated(uint32_t* glo_dist, 
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
            tmp = (((-1-(sqrt((float) d)))/2)*(-1))+0.9999;
            next = (int) tmp;
            i = (cities-2) - (next-1);
            j = (i+2) + (ind-(totIter-((next*(next-1))/2))); 
            ip1 = i+1;
            jp1 = j+1;
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
        num_threads = (num_elems + 1)/2;

        //Reduction on all the local minimum changes found by each thread
        //to find the best minimum change for this climber.
        while(num_threads != num_elems){
            if (idx < num_threads){
                tempRes[idx] = minInd::apply(tempRes[idx],tempRes[idx + num_threads]);
            }
            __syncthreads();

            num_elems = num_threads;
            num_threads = (num_elems + 1)/2;
        }
        ChangeTuple best = minInd::remVolatile(tempRes[0]);
        //Prepare information for swapping
        int temp, swapCities;
        i = best.i + 1;
        j = best.j;
        swapCities = (((j - best.i) + 1)/2) + i; //the ceiling of j/2 plus i
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
    
    int local_opt_cost = sumTourKernel(glo_dist, tour, cities, tempRes);

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
