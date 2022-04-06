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


__device__ int sumTourKernel(uint32_t* glo_dist, unsigned short *lo_tour, int cities, int* result_arr){
    int idx = threadIdx.x;
    int sum = 0;
    int glo_i, glo_ip1; 
    for(int i = idx; i < cities; i += blockDim.x){
        glo_i = lo_tour[i];
        glo_ip1 = lo_tour[i+1];
        sum += glo_dist[glo_i * cities + glo_ip1];
    }
    result_arr[idx] = sum;
    __syncthreads();
    for (int size = blockDim.x /2; size > 0; size /= 2 ){
        if (idx < size) {
            result_arr[idx] += result_arr[idx+size];
        }
        __syncthreads();
    }
    return result_arr[idx];
}

__global__ void twoOptKer(uint32_t* glo_dist, unsigned short *glo_tour, int* glo_is, int* glo_js, int cities, int totIter){
    int block_size = blockDim.x;
    int idx = threadIdx.x;
    int i, j;
    int32_t localMinChange[3];
    extern __shared__ unsigned short totShared[];   //shared memory for both tour, minChange and tempRes
    unsigned short* tour = totShared;               //tour for this climber
    int* tempRes = (int*)&tour[cities+1];           //tempRes holds the best local changes found by each thread
    int* minChange = (int*)&tempRes[3*block_size];  //minChange holds the current best change


    /* Test of shared memory
    int resSize = blockDim.x + cities+1;
    int totSize = resSize+3;
    for(i = threadIdx.x; i < totSize; i += blockDim.x){
        if(i < cities+1){
            tour[i] = glo_tour[i];
            printf("shareTour: %d\n", shaTour[i]);
        }
        else if(i > cities && i < resSize){
            int tmp = (i-(cities+1))*3;
            tempRes[tmp] = float (tmp);
            tempRes[tmp+1] = float (tmp);  
            tempRes[tmp+2] = float (tmp);  
            printf("temp res: fst %f, sec %f, thr %f \n", tempRes[tmp], tempRes[tmp+1], tempRes[tmp+2]);
        }else if(i < totSize){
            int tmpM = (i-resSize)*3;
            minChange[tmpM] = float (tmpM);
            minChange[tmpM+1] = float (tmpM);
            minChange[tmpM+2] = float (tmpM);
            printf("minChange: fst %f, sec %f, thr %f\n", minChange[tmpM], minChange[tmpM+1], minChange[tmpM+2]);
        }
    }*/
    if(idx == 0) {
        printf("check 1 \n");
    }
    //Preparing data for the 2 opt algorithm
    int ip1, jp1, change;
    for(int t = idx; t < cities+2; t += block_size){
        if(t < cities+1){   //initialize tour to shared memory
            tour[t] = glo_tour[t];
        }else{              //initialize minChange to shared memory
            minChange[0] = -1; 
            minChange[1] = 0; 
            minChange[2] = 0;
            //printf("thread num in else: %d \n", t);
        }
    }
    
    __syncthreads();
    if(idx == 0) {
        printf("check 2 \n");
    }
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
            //printf("each threads smallest element: change %d, i %d, j %d \n", localMinChange[0], localMinChange[1], localMinChange[2]);   
        }
        //Write each threads local minimum change (best change found)
        //to the shared array tempRes. 
        if(idx < totIter){
            tempRes[idx*3] = localMinChange[0];
            tempRes[idx*3+1] = localMinChange[1];
            tempRes[idx*3+2] = localMinChange[2];
            //printf("res: change %d, i %d, j %d \n", tempRes[idx*3], tempRes[idx*3+1], tempRes[idx*3+2]);
        }
        __syncthreads();
        if(idx == 0) {
        printf("check 3 \n");
        }
        //Preparation for the reduction on all local minimum changes.
        int num_elems, num_threads;
        if(totIter < block_size){
            num_elems = totIter;
        }else{
            num_elems = block_size;
        }
        num_threads = (num_elems + 1 ) / 2;
        if(idx == 0) {
        printf("check 3 \n");
    }
        //Reduction on all the local minimum changes found by each thread
        //to find the best minimum change for this climber.
        while(1){
            if(idx < num_threads){
                if (idx + num_threads < num_elems){
                    //printf("num_th: %d, num_elem: %d\n",num_threads, num_elems);
                    //printf("threadid: %d, threadid+threadnum: %d\n",idx*3, (idx + num_threads)*3);
                    //printf("first elem: %d, second elem: %d\n",tempRes[idx*3], tempRes[(idx + num_threads)*3]);
                    if (tempRes[idx*3] > tempRes[(idx + num_threads)*3]) {
                        //printf("if statement\n");
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
            if(idx == 0) {
            printf("check 5 \n");
            }
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
        //printf("i: %d, j: %d, swapc: %d\n ", i, j, swapCities);
        //swap
        for(int t = idx + i; t < swapCities; t += block_size){
            //printf("t: %d, swapc: %d\n ", t, swapCities);
            temp = tour[t];
            tour[t] = tour[j - (t - i)];
            tour[j - (t - i)] = temp;
        }
        if(idx < 3){
            minChange[idx] = tempRes[idx];
            //printf("idx: %d, minChange: %d\n ", idx, minChange[idx]);

        }
        __syncthreads();
    }
    int local_opt_cost = sumTourKernel(glo_dist, tour, cities, tempRes);
    printf("idx: %d, local cost: %d\n", idx, tempRes[0]);
    
}

