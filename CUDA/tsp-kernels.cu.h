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

__global__ void twoOptKer(uint32_t *gloDist, unsigned short *gloTour, int cities){
    int i, j, ip1, jp1;
    extern __shared__ unsigned short totShared[]; //shared memory for both tour, minChange and tempResults
    unsigned short* shaTour = totShared;
    float* tempRes = (float*)&totShared[cities+1];
    float* minChange = (float*)&tempRes[3*blockDim.x];
    int resSize = blockDim.x + cities+1;
    int totSize = resSize+3;
    for(i = threadIdx.x; i < ; i += blockDim.x){
        if(i < cities+1){
            shaTour[i] = gloTour[i];
            printf("shareTour: %d\n", shaTour[i]);
        }
        else if(i > cities && i < resSize){
            int tmp = (i-(cities+1))*3;
            tempRes[tmp] = float (tmp);
            tempRes[tmp+1] = float (tmp);  
            tempRes[tmp+2] = float (tmp);  
            printf("temp res: fst %02.f, sec %02.f, thr %02.f \n", tempRes[tmp], tempRes[tmp+1], tempRes[tmp+2]);
        }else if(i < totSize){
            int tmpM = (i-resSize)*3;
            minChange[tmpM] = float (tmpM);
            minChange[tmpM+1] = float (tmpM);
            minChange[tmpM+2] = float (tmpM);
            printf("minChange: fst %02.f, sec %02.f, thr %02.f\n", minChange[tmpM], minChange[tmpM+1], minChange[tmpM+2]);
        }
    }
    __syncthreads();



}