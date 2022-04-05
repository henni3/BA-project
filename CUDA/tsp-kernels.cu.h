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

__global__ void twoOptKer(uint32_t* glo_dist, unsigned short *glo_tour, int* glo_is, int* glo_js, int cities, int totIter){
    int i, j;
    int32_t localMinChange[3];
    extern __shared__ unsigned short totShared[];   //shared memory for both tour, minChange and tempRes
    unsigned short* tour = totShared;               //tour for this climber
    int* tempRes = (int*)&tour[cities+1];           //tempRes holds the best local changes found by each thread
    int* minChange = (int*)&tempRes[3*blockDim.x];  //minChange holds the current best change


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

    //Preparing data for the 2 opt algorithm
    int ip1, jp1, change;
    for(int t = threadIdx.x; t < cities+2; t += blockDim.x){
        if(t < cities+1){   //initialize tour to shared memory
            tour[t] = glo_tour[t];
        }else{              //initialize minChange to shared memory
            minChange[0] = -1; 
            minChange[1] = 0; 
            minChange[2] = 0;
            //printf("thread num in else: %d \n", t);
        }
    }
    localMinChange[0] = 0; 
    localMinChange[1] = 0; 
    localMinChange[2] = 0;
    __syncthreads();

    //Computation for one climber
    while(minChange[0] < 0){
        if(threadIdx.x == 0){
            minChange[0] = 0;
        }
        /***
        The 2 opt move
        Each thread calculates the local changes of the given i and j indexes.
        The i and j index are collected (with a stride of block size) from the 
        global i array and in the global j array to acheive coalesecing.
        ***/
        for(int ind = threadIdx.x; ind < totIter; ind += blockDim.x){
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
        if(threadIdx.x < totIter){
            tempRes[threadIdx.x*3] = localMinChange[0];
            tempRes[threadIdx.x*3+1] = localMinChange[1];
            tempRes[threadIdx.x*3+2] = localMinChange[2];
            //printf("res: change %d, i %d, j %d \n", tempRes[threadIdx.x*3], tempRes[threadIdx.x*3+1], tempRes[threadIdx.x*3+2]);
        }
        __syncthreads();
        
        //Preparation for the reduction on all local minimum changes.
        int num_elems, num_threads;
        if(totIter < blockDim.x){
            num_elems = totIter;
        }else{
            num_elems = blockDim.x;
        }
        num_threads = (num_elems + 1 ) / 2;

        //Reduction on all the local minimum changes found by each thread
        //to find the best minimum change for this climber.
        while(threadIdx.x < num_threads){
            if (threadIdx.x + num_threads < num_elems){
                //printf("num_th: %d, num_elem: %d\n",num_threads, num_elems);
                //printf("threadid: %d, threadid+threadnum: %d\n",threadIdx.x*3, (threadIdx.x + num_threads)*3);
                //printf("first elem: %d, second elem: %d\n",tempRes[threadIdx.x*3], tempRes[(threadIdx.x + num_threads)*3]);
                if (tempRes[threadIdx.x*3] > tempRes[(threadIdx.x + num_threads)*3]) {
                    //printf("if statement\n");
                    tempRes[threadIdx.x*3] = tempRes[(threadIdx.x + num_threads)*3];
                    tempRes[threadIdx.x*3 + 1] = tempRes[(threadIdx.x + num_threads)*3 + 1];
                    tempRes[threadIdx.x*3 + 2] = tempRes[(threadIdx.x + num_threads)*3 + 2];
                }
                else if (tempRes[threadIdx.x*3] == tempRes[(threadIdx.x + num_threads)*3]){
                    if (tempRes[threadIdx.x*3 + 1] > tempRes[(threadIdx.x + num_threads)*3 +1 ]){
                        tempRes[threadIdx.x*3] = tempRes[(threadIdx.x + num_threads)*3 ];
                        tempRes[threadIdx.x*3 + 1] = tempRes[(threadIdx.x + num_threads)*3 + 1];
                        tempRes[threadIdx.x*3 + 2] = tempRes[(threadIdx.x + num_threads)*3 + 2];
                    }
                    else if (tempRes[threadIdx.x*3 +1] == tempRes[(threadIdx.x + num_threads)*3 + 1]){
                        if (tempRes[threadIdx.x*3 + 2] > tempRes[(threadIdx.x + num_threads)*3] +2){
                            tempRes[threadIdx.x*3] = tempRes[(threadIdx.x + num_threads)*3 ];
                            tempRes[threadIdx.x*3 + 1] = tempRes[(threadIdx.x + num_threads)*3 + 1];
                            tempRes[threadIdx.x*3 + 2] = tempRes[(threadIdx.x + num_threads)*3 + 2];
                        }
                    }
                }
            }
            //printf("thread id: %d\n", threadIdx.x);
            __syncthreads();

            num_elems = num_threads;
            num_threads= (num_elems + 1)/ 2;
            if(num_threads == num_elems){
                break;
            }
        }
        printf("hej\n");
        //Prepare information for swapping
        int temp, swapCities;
        i = tempRes[1] + 1;
        j = tempRes[2];
        swapCities = (((tempRes[1] - tempRes[2]) + 1) / 2) + i; //the ceiling of j/2 plus i
        //swap
        for(int t = threadIdx.x + i; t < swapCities; t += blockDim.x){
            printf("t: %d, swapc: %d\n ", t, swapCities);
            temp = tour[t];
            tour[t] = tour[j - (t - i)];
            tour[j - (t - i)] = temp;
        }
        __syncthreads();
        if(threadIdx.x < (cities+1)){
            printf("tID: %d, tourElem: %d\n", threadIdx.x, tour[threadIdx.x]);
        }
        
        

    }
}