__global__ void
replicate0(int tot_size, char* flags_d) {
    int glb_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(glb_id < tot_size){
        flags_d[glb_id]=0;
    }
}

__global__ void
mkFlags(int mat_rows, int* mat_shp_sc_d, char* flags_d) {
    int glb_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(glb_id < mat_rows){
        if(glb_id == 0){
             flags_d[0]=1;
        }else{
            flags_d[mat_shp_sc_d[glb_id-1]]=1;
        }
        
    }
}