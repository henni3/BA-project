#ifndef TRANSPOSE_KERS
#define TRANSPOSE_KERS

typedef unsigned int uint32_t; 

/**
 * Naive memcpy kernel, for the purpose of comparing with
 * a more "realistic" bandwidth number.
 */
__global__ void naiveMemcpy(float* d_out, float* d_inp, const uint32_t N) {
    uint32_t gid = blockIdx.x * blockDim.x + threadIdx.x;
    if(gid < N) {
        d_out[gid] = d_inp[gid];
    }
}


// widthA = heightB
template <class T> 
__global__ void matTransposeKer(T* A, T* B, int heightA, int widthA) {

  int gidx = blockIdx.x*blockDim.x + threadIdx.x;
  int gidy = blockIdx.y*blockDim.y + threadIdx.y; 

  if( (gidx >= widthA) || (gidy >= heightA) ) return;

  B[gidx*heightA+gidy] = A[gidy*widthA + gidx];
}

// blockDim.y = T; blockDim.x = T
// each block transposes a square T
template <class ElTp, int T> 
__global__ void matTransposeTiledKer(ElTp* A, ElTp* B, int heightA, int widthA) {
  extern __shared__ char sh_mem1[];
  volatile ElTp* tile = (volatile ElTp*)sh_mem1;
  //__shared__ float tile[T][T+1];

  int x = blockIdx.x * T + threadIdx.x;
  int y = blockIdx.y * T + threadIdx.y;

  if( x < widthA && y < heightA )
      tile[threadIdx.y*(T+1) + threadIdx.x] = A[y*widthA + x];

  __syncthreads();

  x = blockIdx.y * T + threadIdx.x; 
  y = blockIdx.x * T + threadIdx.y;

  if( x < heightA && y < widthA )
      B[y*heightA + x] = tile[threadIdx.x*(T+1) + threadIdx.y];
}


__global__ void 
origProg(float* A, float* B, unsigned int N) {
    unsigned int gid = (blockIdx.x * blockDim.x + threadIdx.x);
    if(gid >= N) return;
    float accum = 0.0;
    unsigned int thd_offs = gid * 64;

    for(int j=0; j<64; j++) {
        float tmpA  = A[thd_offs + j];
        accum = sqrt(accum) + tmpA*tmpA;
        B[thd_offs + j]  = accum;
    }
}

/**
 * To achieve coalesced access we read from the transposed of A, named Atr,
 * and write the result in the transposed of B, named Btr. The total number
 * of threads is N, i.e., the number of rows of A and B (and the number of
 * columns of Atr and Btr).
 *
 * Task3.a: fill in the implementation of this kernel,
 *          such that all accesses are coalesced.
 */
__global__ void 
transfProg(float* Atr, float* Btr, unsigned int N) {
    unsigned int gid = (blockIdx.x * blockDim.x + threadIdx.x);
    if(gid >= N) return;
    float accum = 0.0;
    for(int j=0; j<64; j++) {
        float tmpA  = Atr[j*N+gid];
        accum = sqrt(accum) + tmpA*tmpA;
        Btr[j*N+gid] = accum;
    }
}

/**
 * Task 3.b: fill in the implementation of this function, which is
 *           used to copy "CHUNK*blockDim.x" elements between shared
 *           and global memory in coalesced fashion---please see how
 *           it is called in kernel "optimProg". 
 *           This function is supposed to compute:
 *           "*loc_ind": the shared-memory index where the current 
 *                       thread should write/read and
 *           "*glb_ind": the global-memory index where the current
 *                       thread should read/write (from/to arrays
 *                       A and B)
 * 
 * Input arguments are:
 * "block_offs = blockIdx.x * blockDim.x * 64" is the offset of the
 *              current block
 * "jj" is the count of the virtualization loop "jj=0; jj<64; jj+=CHUNK"
 * "j"  takes values "j=0; j<CHUNK; j++"
 * "chunk_lane = threadIdx.x % CHUNK"
 * "chunk_id   = threadIdx.x / CHUNK"
 * The assumption is that the block size is a multiple of CHUNK---this
 *   is verified in the host code.
 *
 * HINTS:
 * "*loc_ind" is relatively easy: you can see the shared memory as a two
 *            dimensional array of outer size [CHUNK] and inner size 
 *            [blockDim.x], where threadIdx.x spans the inner dimension.
 *            What can the outer-dimension index be?
 *
 * "*glb_ind" is more complicated because you have to look at it as a 2D
 *            array of inner size 64. The idea is that groups of CHUNK
 *            consecutive threads read from the same row (so that the
 *            accesses are coalesced). Please reason for each component
 *            individually:
 *    1. you probably need to add the block_offs
 *    2. If a group of CHUNK threads accesses the same row, then with each
 *       new "j", how many rows do you need to jump over?
 *    3. "jj" and "chunk_lane" probably contribute to the inner dimension
 *    4. how about "chunk_id"? Does it contribute to the inner dimension
 *       or to the outer dimension? 
 *       (i.e., do you need to multiply it with 64 or not?)
 *            
 */ 
template<int CHUNK>
__device__ inline void 
glb2shmem( uint32_t block_offs, uint32_t jj, uint32_t j
         , uint32_t chunk_lane, uint32_t chunk_id // input
         , uint32_t* glb_ind, uint32_t* loc_ind // the result  
) {
    *loc_ind = 0;
    *glb_ind = 0;
}


/**
 * A and B are the original matrices having N rows and 64 elements per row.
 * This kernel is executed by N threads, each thread processing independently
 *   its own row. However, thread cooperate in shared memory---they use the
 *   shared memory (of length CHUNK*blockDim.x) as a staging buffer to optimize
 *   global-memory coalescing.
 */
template<int CHUNK>
__global__ void 
optimProg(float* A, float* B, unsigned int N) {
    extern __shared__ float sh_mem[]; // length: CHUNK * blockDim.x
    
    uint32_t block_offs = blockIdx.x * blockDim.x * 64;
    uint32_t gid = blockIdx.x * blockDim.x + threadIdx.x;
    float accum = 0.0;

    uint32_t chunk_lane = threadIdx.x % CHUNK;
    uint32_t chunk_id   = threadIdx.x / CHUNK;

    for(int jj=0; jj<64; jj+=CHUNK) {
        uint32_t count_j = min(jj+CHUNK, 64);
       
        // copy CHUNK elements per thread from global to shared
        // memory in coalesced fashion (i.e., coalescing degree is CHUNK)
        for(int j=0; j<CHUNK; j++) {
            uint32_t loc_ind, glb_ind;
            glb2shmem<CHUNK>( block_offs, jj, j, chunk_lane, chunk_id, &glb_ind, &loc_ind );
            if(glb_ind < N*64) {
                sh_mem[loc_ind] = A[glb_ind];
            } 
        }
        __syncthreads();

        if(gid < N)
        for(int j=jj; j<count_j; j++) {
            // look at shared memory as a float[blockDim.x][CHUNK] array
            // each thread processes its CHUNK consecutive elements from
            // shared memory and updates it in-place with the result
            uint32_t loc_ind = threadIdx.x * CHUNK + j-jj;
            float tmpA = sh_mem[loc_ind];
            accum = sqrt(accum) + tmpA*tmpA;
            sh_mem[loc_ind] = accum;
        }
        __syncthreads();

        // copy CHUNK elements per thread from shared to global
        // memory in coalesced fashion (i.e., coalescing degree is CHUNK)
        for(int j=0; j<CHUNK; j++) {
            uint32_t loc_ind, glb_ind;
            glb2shmem<CHUNK>( block_offs, jj, j, chunk_lane, chunk_id, &glb_ind, &loc_ind );
            if(glb_ind < N*64) {
                B[glb_ind] = sh_mem[loc_ind];
            } 
        }
        __syncthreads();
    }
}

#endif
