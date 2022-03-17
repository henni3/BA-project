
__global__ void mkIndShp(int* index_shp_d, int len){
    int glb_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(glb_id < len){
        index_shp_d[glb_id] = len - glb_id;
    }
}

__global__ void replicate0(int totIter, char* flags_d) {
    int glb_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(glb_id < totIter){
        flags_d[glb_id]=0;
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