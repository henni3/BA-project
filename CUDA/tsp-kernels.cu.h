#include "tsp-ker-helper.cu.h"
#include <sys/time.h>
#include "constants.cu.h"
#include <stdio.h>
/* 
 * This file is kernel functions from where the 
 * version 1, 2, 3 and 4 are executed.
 */

__global__ void twoOptKer(uint32_t* glo_dist, 
                          unsigned short *glo_tours, 
                          int* glo_is,
                          int* glo_result, 
                          int cities, 
                          int totIter){
    int block_size = blockDim.x;
    int idx = threadIdx.x;
    int i, j, change, ip1, jp1;
    ChangeTuple localMinChange;
    ChangeTuple maxValue = ChangeTuple(INT_MAX, USHRT_MAX, USHRT_MAX);
    
    //shared memory for both tour, minChange and tempRes
    extern __shared__ unsigned char totShared[];             
    volatile ChangeTuple* tempRes = (volatile ChangeTuple*)&totShared;       //tempRes holds the best local changes found by each thread
    volatile ChangeTuple* minChange = tempRes + block_size;                 //minChange holds the current best change
    volatile unsigned short* tour =
                 (volatile unsigned short*)(minChange + 1);                 //tour for this climber

    if(minChange == NULL){
        printf("pointer error\n");
    }

    //copy global tour to shared memory
    for(int t = idx; t < cities+1; t += block_size){
        tour[t] = glo_tours[blockIdx.x * (cities+1) + t];
    }
    
    //initialize minChange to shared memory
    if(idx == 0){
        minChange[0] = ChangeTuple();
        minChange[0].change = -1;
    }
    
    __syncthreads();
    //Computation for one climber
    while(minChange[0].change < 0){
        __syncthreads();
        if(idx == 0){
            minChange[0] = ChangeTuple();
        }
        //reset each threads local min change
        localMinChange = ChangeTuple();
        
        /***
        The 2 opt move
        Each thread calculates the local changes of the given i and j indexes.
        The i and j index are collected (with a stride of block size) from the 
        global i array and in the global j array to acheive coalesecing.
        ***/
        for(int ind = idx; ind < totIter; ind += block_size){
            int num = glo_is[ind]; // 4
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
        //to the shared array tempRes. If there are threads does not hold a
        //calculated local minimum change (i.e. totIter < block_size), then
        //write the maximum values to tempRes.
        if(idx < totIter){
            tempRes[idx] = ChangeTuple(localMinChange);
        }else{
            tempRes[idx] = ChangeTuple(maxValue);
        }
        __syncthreads();
        
        //Reduction on all the local minimum changes found by each thread
        //to find the best minimum change for this climber.
        reduceLocalMinChange(block_size, tempRes);

        //The best local minimum changes found is stored in variable best.
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
        if(idx == 0){
            minChange[0].change = tempRes[0].change;
            minChange[0].j = tempRes[0].j;
            minChange[0].i = tempRes[0].i;
        }
        __syncthreads();
    }

    
    int local_opt_cost = sumTourKernel(glo_dist, tour, cities, tempRes);


    //copy best local tour from shared memory to global memory
    for(int t = idx; t < cities+1; t += block_size){
        glo_tours[blockIdx.x * (cities+1) + t] = tour[t];
    }

    
    //Writing local optimum results to global memory
    if(idx == 0){
        glo_result[blockIdx.x * 2] = local_opt_cost;
        glo_result[blockIdx.x * 2+1] = blockIdx.x;
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
    ChangeTuple maxValue = ChangeTuple(INT_MAX, USHRT_MAX, USHRT_MAX);

    extern __shared__ unsigned char totShared[];             //shared memory for both tour, minChange and tempRes
    volatile ChangeTuple* tempRes = (volatile ChangeTuple*)&totShared;       //tempRes holds the best local changes found by each thread
    volatile ChangeTuple* minChange = tempRes + block_size; //minChange holds the current best change (cos?)
    volatile uint32_t* shared_Dist = (volatile uint32_t*) (minChange + 1);      
    volatile unsigned short* tour =
                 (volatile unsigned short*)(shared_Dist + cities * cities);  //tour for this climber
    if(minChange == NULL){
        printf("pointer error\n");
    }

    ChangeTuple tmpMinChange = ChangeTuple();

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
            tempRes[idx] = ChangeTuple(localMinChange);
        }else{ //else pad with neutral element
            tempRes[idx] = ChangeTuple(maxValue);
        }
        __syncthreads();
        
        //Reduction on all the local minimum changes found by each thread
        //to find the best minimum change for this climber.
        reduceLocalMinChange(block_size, tempRes);
        
        //The best local minimum changes found is stored in variable best.
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
        
        tmpMinChange.change = tempRes[0].change;
        tmpMinChange.j = tempRes[0].j;
        tmpMinChange.i = tempRes[0].i;
        __syncthreads();
    } 
    int local_opt_cost = sumTourKernel100(shared_Dist, tour, cities, tempRes);

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
    int change, d, i, j, ip1, jp1, block_size, idx;
    unsigned int curr;
    block_size = blockDim.x;
    idx = threadIdx.x;
    ChangeTuple localMinChange;
    ChangeTuple maxValue = ChangeTuple(INT_MAX, USHRT_MAX, USHRT_MAX);

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
    
    //initialize minChange to shared memory
    if(idx == 0){
        minChange[0] = ChangeTuple();
        minChange[0].change = -1;
    }
    __syncthreads();

    //Computation for one climber
    while(minChange[0].change < 0){
        __syncthreads();
        if(idx == 0){
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
            tmp = ((-1+(sqrt((float) d)))/2)+0.9999;
            curr = (int) tmp;
            i = (cities-2) - curr;
            j = (i+2) + (ind-(totIter-((curr*(curr+1))/2))); 
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
            tempRes[idx] = ChangeTuple(localMinChange);
        }else{ //else pad with neutral element
            tempRes[idx] = ChangeTuple(maxValue);
        }
        __syncthreads();
        
        //Reduction on all the local minimum changes found by each thread
        //to find the best minimum change for this climber.
        reduceLocalMinChange(block_size, tempRes);
        
        ChangeTuple best = minInd::remVolatile(tempRes[0]);
        //Prepare information for swapping
        int temp, swapCities;
        i = best.i + 1;
        j = best.j;
        swapCities = (((j - best.i) + 1)>>1) + i; //the ceiling of j/2 plus i
        //swap
        for(int t = idx + i; t < swapCities; t += block_size){
            temp = tour[t];
            tour[t] = tour[j - (t - i)];
            tour[j - (t - i)] = temp;
        }
        if(idx == 0){
            minChange[0].change = tempRes[0].change;
            minChange[0].j = tempRes[0].j;
            minChange[0].i = tempRes[0].i;
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

__global__ void twoOptKer_test(uint32_t* glo_dist, 
                          unsigned short *glo_tours,
                          int* glo_result,
                          int* counter, 
                          int cities, 
                          int totIter){
    int block_size = blockDim.x;
    int idx = threadIdx.x;
    int change, d, i, j, ip1, jp1;
    unsigned int curr;
    ChangeTuple localMinChange;
    ChangeTuple maxValue = ChangeTuple(INT_MAX, USHRT_MAX, USHRT_MAX);
    
    //shared memory for both tour, minChange and tempRes
    extern __shared__ unsigned char totShared[];             
    volatile ChangeTuple* tempRes = (volatile ChangeTuple*)&totShared;       //tempRes holds the best local changes found by each thread
    volatile ChangeTuple* minChange = tempRes + block_size;                 //minChange holds the current best change
    volatile int* while_block = (volatile int*) (minChange + 1);
    volatile unsigned short* tour =
                 (volatile unsigned short*)(while_block + block_size);                 //tour for this climber
    //printf("test1.12 \n");

    if(minChange == NULL){
        printf("pointer error\n");
    }

    //extern__shared__ int while_block[1024]; //max blocksize, might be wasted for small inputs

    // Init of counter array 
    while_block[idx] = 0;
    

    //copy global tour to shared memory
    for(int t = idx; t < cities+1; t += block_size){
        tour[t] = glo_tours[blockIdx.x * (cities+1) + t];
        // cities + 1 * 2 bytes
    }
    
    //initialize minChange to shared memory
    if(idx == 0){
        minChange[0] = ChangeTuple();
        minChange[0].change = -1;
    }
    
    __syncthreads();
    //Computation for one climber
    while(minChange[0].change < 0){
        while_block[idx]++;
        __syncthreads();
        if(idx == 0){
           //repeats++;
            minChange[0] = ChangeTuple();
        }
        //reset each threads local min change
        localMinChange = ChangeTuple();
        
        /***
        The 2 opt move
        Each thread calculates the local changes of the given i and j indexes.
        The i and j index are collected (with a stride of block size) from the 
        global i array and in the global j array to acheive coalesecing.
        ***/
        float tmp;
        for(int ind = idx; ind < totIter; ind += block_size){
            d = 1-(4*(-2*(totIter-ind)));
            tmp = ((-1+(sqrt((float) d)))/2)+0.9999;
            curr = (int) tmp;
            i = (cities-2) - curr;
            j = (i+2) + (ind-(totIter-((curr*(curr+1))/2))); 
            ip1 = i+1;
            jp1 = j+1;
            change = glo_dist[tour[i]*cities+tour[j]] + 
                    glo_dist[tour[ip1]*cities+tour[jp1]] -
                    (glo_dist[tour[i]*cities+tour[ip1]] +
                    glo_dist[tour[j]*cities+tour[jp1]]);
            // 4 * 4
            //Each thread shall hold the best local change found
            ChangeTuple check = ChangeTuple(change,(unsigned short)i, (unsigned short) j);
            localMinChange = minInd::apply(localMinChange,check);
        }
         // for loop = 4*4 * totiter bytes
        //Write each threads local minimum change (best change found)
        //to the shared array tempRes. If there are threads does not hold a
        //calculated local minimum change (i.e. totIter < block_size), then
        //write the maximum values to tempRes.
        if(idx < totIter){
            tempRes[idx] = ChangeTuple(localMinChange);
        }else{
            tempRes[idx] = ChangeTuple(maxValue);
        }
        __syncthreads();
        
        //Reduction on all the local minimum changes found by each thread
        //to find the best minimum change for this climber.
        reduceLocalMinChange(block_size, tempRes);

        //The best local minimum changes found is stored in variable best.
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
        if(idx == 0){
            minChange[0].change = tempRes[0].change;
            minChange[0].j = tempRes[0].j;
            minChange[0].i = tempRes[0].i;
        }
        __syncthreads();
    }
    __syncthreads();

     // reduceLocalCounter(block_size, while_block);
     // while loop = 4*4 * totiter * While_iters
    
    int local_opt_cost = sumTourKernel(glo_dist, tour, cities, tempRes); //cosmin do we do this multiple times ( uneccessary computation = threads -1 )

    // 4 * cities 


    //copy best local tour from shared memory to global memory
    for(int t = idx; t < cities+1; t += block_size){
        glo_tours[blockIdx.x * (cities+1) + t] = tour[t];
    }
    // 4 * cities+1 bytes
    
    //Writing local optimum results to global memory
    if(idx == 0){
        glo_result[blockIdx.x * 2] = local_opt_cost;
        glo_result[blockIdx.x * 2+1] = blockIdx.x;
        counter[blockIdx.x] = while_block[0];
        //counter[blockIdx.x * 2 + 1]  = blockIdx.x;
        //printf("number of while iters in block %d is : %d  \n", blockIdx.x, while_block[0] );
    }
    // 4 * 2 

    //total work for one block is 6 * cities + 14 + (16 * totiter * while)
}
