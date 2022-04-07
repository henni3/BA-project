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

//Random tour generator basen on SPLASH-2 code
/*__global__ void createTour (unsigned short* iniTour, int cities, int tourOffset){
    int currTour, rand, mult, add, mask, idx, to, temp;
    idx = threadIdx.x;
    currTour = idx + blockIdx.x * blockDim.x;
    rand = currTour + tourOffset;
    mult = 1103515245;
    add = 12345;
    mask = 0x7fffffff;
    for(int t = idx; t < cities; t += block_size){
        rand = (mult * rand + add) & mask;
        to = rand % cities;
        if (to <= 0){
            to = 1;
        }
        temp = tour[t];
        tour[t] = tour[to];
        tour[to] = temp;
    }

} */


__global__ void twoOptKer(uint32_t* glo_dist, 
                          unsigned short *glo_tour, 
                          int* glo_is, int* glo_js, 
                          int cities, 
                          int totIter){
    int block_size = blockDim.x;
    int idx = threadIdx.x;
    int i, j;
    int32_t localMinChange[3];
    //printf("hej1 \n");
    extern __shared__ unsigned short totShared[];   //shared memory for both tour, minChange and tempRes
    unsigned short* tour = totShared;               //tour for this climber
    //printf("hej2 \n");
    int* tempRes = (int*)&tour[cities+1];           //tempRes holds the best local changes found by each thread
    //printf("hej3 \n");
    int* minChange = (int*)&tempRes[3*block_size];  //minChange holds the current best change
    if(minChange == NULL){
        printf("pointer error\n");
    }


    //Test of shared memory
    int resSize = blockDim.x + cities+1;
    int totSize = resSize+3;
    for(i = threadIdx.x; i < totSize; i += blockDim.x){
        if(i < cities+1){
            tour[i] = glo_tour[i];
            printf("shareTour: %d\n", tour[i]);
        }else if(i > cities && i < resSize){
            printf("i: %d", i);
            int tmp = (i-(cities+1))*3;
            printf("tmp: %d", tmp);
            tempRes[tmp] = tmp;
            tempRes[tmp+1] = tmp;  
            tempRes[tmp+2] = tmp;  
            printf("temp res: fst %d, sec %d, thr %d \n", tempRes[tmp], tempRes[tmp+1], tempRes[tmp+2]);
        }else if(i < totSize){
            int tmpM = (i-resSize)*3;
            minChange[tmpM] = tmpM;
            minChange[tmpM+1] = tmpM;
            minChange[tmpM+2] = tmpM;
            printf("minChange: fst %d, sec %d, thr %d\n", minChange[tmpM], minChange[tmpM+1], minChange[tmpM+2]);
        }
    }
    if(idx == 0) {
        printf("check 1 \n");
    }
    //Preparing data for the 2 opt algorithm
    int ip1, jp1, change;
    //initialize tour to shared memory
    //printf("hej1 \n");
    for(int t = idx; t < cities+1; t += block_size){
        tour[t] = glo_tour[t];
        printf("idx %d, tour: %d\n", t, tour[t]);
    }
    if(idx == 0) {
        printf("before Change\n");
    }
    minChange[0] = 0;
    if(idx == 0) {
        printf("before if, idx %d \n ", idx);

    }
    if(idx == 0){
        printf("in if, with thread id %d = 0, with minchange %d \n ", idx, minChange[0]);
        //initialize minChange to shared memory
        minChange[0] = -1; 
        //printf("gets here 1 \n");
        minChange[1] = 0; 
        //printf("gets here 2 \n");
        minChange[2] = 0;
        printf("after minchanges, where the values are min_1 %d, min2 %d, min3 %d \n",minChange[0],minChange[1],minChange[2]);
    }
    //printf("before sync \n ");
    
    if(idx == 0) {
        printf("check 2 \n");
    }
    __syncthreads();
    //Computation for one climber
    if(idx == 0) {
        printf("before while \n");
    }
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
            if(ind == 0) {
                printf("in for\n");
            }
            i = glo_is[ind];
            j = glo_js[ind] + i + 2;
            printf("tour in for. i: %d. j: %d\n", i, j);

            if(ind == 0) {
                printf("forbi j\n");
            }
            ip1 = i+1;
            if(ind == 0) {
                printf("forbi ip1\n");
            }
            jp1 = j+1;
            if(ind == 0) {
                printf("forbi jp1\n");
            }
            change = glo_dist[tour[i]*cities+tour[j]] + 
                    glo_dist[tour[ip1]*cities+tour[jp1]] -
                    (glo_dist[tour[i]*cities+tour[ip1]] +
                    glo_dist[tour[j]*cities+tour[jp1]]);
            if(ind == 0) {
                printf("after change \n");
           }
            //Each thread shall hold the best local change found
            if(change < localMinChange[0]){  
                if(ind == 0) {
                    printf("change local \n");
                }
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

