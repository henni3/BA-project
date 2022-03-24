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

__global__ void twoOptKer(uint32_t *gloDist, uint16_t *gloTour, int cities){
    int i, j;
    extern __shared__ uint16_t sharedTourDist[]; //shared memory for both cities and tour
    uint16_t* shaTour = sharedTourDist;
    uint32_t* shaDist = (uint32_t*)&shaTour[cities+1];
    for(i = threadIdx.x; i < cities * cities; i += blockDim.x){
        if(i < cities+1){
            shaTour[i] = gloTour[i];
            shaDist[i] = gloDist[i];
            printf("shareTour: %d\n", shaTour[i]);
        }else{
            shaDist[i] = gloDist[i];
        }
        printf("shareDist: %d\n", shaDist[i]);
    }
    __syncthreads();



}