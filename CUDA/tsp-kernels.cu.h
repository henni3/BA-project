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
    int i, j, ip1, jp1, change;
    int32_t localMinChange[3];
    extern __shared__ unsigned short totShared[]; //shared memory for both tour, minChange and tempRes
    unsigned short* tour = totShared;
    int* tempRes = (int*)&tour[cities+1];
    int* minChange = (int*)&tempRes[3*blockDim.x];


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

    for(int t = threadIdx.x; t < cities+2; t += blockDim.x){
        if(t < cities+1){   //initialize tour to shared memory
            tour[t] = glo_tour[t];
        }else{              //initialize minChange to shared memory
            minChange[0] = -1; 
            minChange[1] = 0; 
            minChange[2] = 0;
            printf("thread num in else: %d \n", t);
        }
    }
    localMinChange[0] = 0; 
    localMinChange[1] = 0; 
    localMinChange[2] = 0;
    __syncthreads();

    while(minChange[0] < 0){
        if(threadIdx.x == 0){
            minChange[0] = 0;
        }
        for(int ind = threadIdx.x; ind < totIter; ind += blockDim.x){
            i = glo_is[ind];
            j = glo_js[ind] + i + 2;
            ip1 = i+1;
            jp1 = j+1;
            change = glo_dist[tour[i]*cities+tour[j]] + 
                    glo_dist[tour[ip1]*cities+tour[jp1]] -
                    (glo_dist[tour[i]*cities+tour[ip1]] +
                    glo_dist[tour[j]*cities+tour[jp1]]);
            if(change < localMinChange[0]){ 
                    localMinChange[0] = change; 
                    localMinChange[1] = i; 
                    localMinChange[2] = j;
            }
            printf("each threads smallest element: change %d, i %d, j %d \n", localMinChange[0], localMinChange[1], localMinChange[2]);   
        }
        for(int t = threadIdx.x; t < totIter; t += blockDim.x){
            int a = tempRes[t*3] = localMinChange[0];
            int b = tempRes[t*3+1] = localMinChange[1];
            int c = tempRes[t*3+2] = localMinChange[2];
            printf("res: change %d, i %d, j %d \n", a, b, c);
        }
    }
        __syncthreads();

}