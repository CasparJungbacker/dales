#include <iostream>
#include <cuda.h>

/* TODO: make this dynamic */
static constexpr int BLOCK_DIM { 512 };

template<typename T>
__global__ void reduction_scalar( T* input, T* output, int n) {
    int idx_g = blockIdx.x * (blockDim.x * 2) + threadIdx.x;
    int idx_l = threadIdx.x;

    __shared__ T smem[BLOCK_DIM];

    if (idx_g + blockDim.x < n) {
       smem[idx_l] = input[idx_g] + input[idx_g + blockDim.x]; 
    } else {
       smem[idx_l] = static_cast<T>(0);
    }

    __syncthreads();

    for (int stride = BLOCK_DIM / 2; stride >= 1; stride >>= 1) {
	    if (idx_l < stride && idx_g + stride < n) {
                smem[idx_l] += smem[idx_l + stride]; 
	}
	__syncthreads();
    }

    if (idx_l == 0) {
        atomicAdd(output, smem[0]);
    }
}

template<typename T>
__global__ void reduction_interior(T* input, T* output, int nh, int itot, int jtot) {

    int tid = threadIdx.x;

    // Compute i and j indices
    int i = nh + threadIdx.x;
    int j = nh + blockIdx.y;

    int gid = i + j * (itot + 2 * nh);

    extern __shared__ T smem[];

    if (i <= itot + nh && j <= jtot + nh) {
        smem[tid] = input[gid];
    } else {
        smem[tid] = static_cast<T>(0);
    }

    __syncthreads();

    for (int stride = itot / 2; stride >= 1; stride >>= 1) {
        if (tid < stride) {
            smem[tid] += smem[tid + stride];
        }
        __syncthreads();
    }

    if (tid == 0) {
        atomicAdd(output, smem[0]);
    }

}

// Interface part, which we link to the Fortran code
extern "C" {
    
    void reduction_scalar_float(float* input, float* output, int n, cudaStream_t stream) {
        int num_blocks = (n + BLOCK_DIM - 1) / (2 * BLOCK_DIM);
        reduction_scalar<float><<<num_blocks, BLOCK_DIM, 0, stream>>>(input, output, n); 
    }

    void reduction_scalar_double(double* input, double* output, int n, cudaStream_t stream) {
        int num_blocks = (n + BLOCK_DIM - 1) / (2 * BLOCK_DIM);  
        reduction_scalar<double><<<num_blocks, BLOCK_DIM, 0, stream>>>(input, output, n); 
    }

    void reduction_2d_float(float* input, float* output, int nh, int itot, int jtot) {
        dim3 grid(1, jtot, 1);
        dim3 block(itot, 1, 1);
	size_t smem_size = itot * sizeof(float);
        reduction_interior<float><<<grid, block, smem_size>>>(input, output, nh, itot, jtot);
        cudaDeviceSynchronize();
    }

    /*
    void reduction_2d_double(double* input, double* output, int nh, int itot, int jtot) {
        dim3 grid(1, jtot + 2 * nh);
        dim3 block(itot + 2 * nh, 1);
        reduction_interior<double><<<grid, block>>>(input, output, nh, itot, jtot);
        cudaDeviceSynchronize();
    }
    */
}

int main() {

    const int itot = 512;
    const int jtot = 512;
    const int nh = 1;

    const int size = (itot + 2*nh)*(jtot + 2*nh);
    const int bytes = size * sizeof(float);

    float* host_input = new float[(itot + 2*nh)*(jtot + 2*nh)];
    float* host_output = new float;

    for (int j = 0; j < jtot + 2 * nh; j++) {
        for (int i = 0; i < itot + 2 * nh; i++) {
            if (i >= nh && i < itot + nh && j >= nh && j < jtot + nh) {
                host_input[i + j*(itot+2*nh)] = static_cast<float>(1);
            } else {
                host_input[i + j*(itot+2*nh)] = static_cast<float>(-5000);
            }
        }
    }

    float* dev_input;
    float* dev_output;

    cudaMalloc(&dev_input, bytes);
    cudaMalloc(&dev_output, sizeof(float));

    cudaMemcpy(dev_input, host_input, bytes, cudaMemcpyHostToDevice);

    reduction_2d_float(dev_input, dev_output, nh, itot, jtot);

    cudaMemcpy(host_output, dev_output, sizeof(float), cudaMemcpyDeviceToHost);

    std::cout << "Result is " << *host_output << std::endl;

    delete[] host_input;
    delete host_output;
    cudaFree(dev_input);
    cudaFree(dev_output);

    return 0;
}
