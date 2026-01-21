#include <iostream>
#include <cuda.h>

inline int next_pow_of_2(int x) { return pow(2, ceil(log(x) / log(2))); };

template<typename T>
__global__ void reduction_scalar(T* input, T* output, int n) {
    int idx_l = threadIdx.x;
    int idx_g = blockIdx.x * (blockDim.x * 2) + threadIdx.x;

    extern __shared__ int smem_pointer[];
		T* smem = reinterpret_cast<T*>(smem_pointer);

    if (idx_g + blockDim.x < n) {
        smem[idx_l] = input[idx_g] + input[idx_g + blockDim.x]; 
    } else {
				smem[idx_l] = static_cast<T>(0);
    }

    __syncthreads();

    for (int stride = blockDim.x / 2; stride >= 1; stride >>= 1) {
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
__global__ void reduction_scalar_halo(T* input, T* output, int nh, int itot, int jtot) {
    int i = threadIdx.x;
    int j = blockIdx.y;
    int ij  = i + j * (itot + 2 * nh) + nh * (1 + itot + 2 * nh);

		extern __shared__ int smem_pointer[];
		T* smem = reinterpret_cast<T*>(smem_pointer);

    if (i < itot / 2) {
        smem[i] = input[ij] + input[ij + itot / 2];
    } else {
        smem[i] = static_cast<T>(0);
    }

    __syncthreads();

    for (int stride = blockDim.x; stride >= 1; stride >>= 1) {
        if (i < stride) {
            smem[i] += smem[i + stride];
        }
        __syncthreads();
    }

    if (i == 0) {
        atomicAdd(output, smem[0]);
    }

}

template<typename T>
__global__ void reduction_profile_halo(T* input, T* output, int nh, int itot, int jtot, int ktot) {
		int i = threadIdx.x;
		int j = blockIdx.y;
		int k = blockIdx.z;
		size_t ijk = i + j * (itot + 2 * nh) + nh * (1 + itot + 2 * nh) + k * ((itot + 2 * nh) * (jtot + 2 * nh));
		
		extern __shared__ int smem_pointer[];
		T* smem = reinterpret_cast<T*>(smem_pointer);

		T sum = static_cast<T>(0);

		for (unsigned int tile = 0; tile < itot / blockDim.x + 1; tile++) {
				unsigned int offset = tile * blockDim.x;
				if (i + offset < itot) {
						sum += input[ijk + offset]; 
				}	
		}

		smem[i] = sum;
		
		__syncthreads();

		// TODO: unroll and use warp primitives
		for (unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
				if (i < stride) {
						smem[i] += smem[i + stride];
				}
				__syncthreads();
		}

		if (i == 0) {
				atomicAdd(&output[k], smem[0]);
		}
}

// Interface part, which we link to the Fortran code
extern "C" {
		
    void reduction_scalar_float(float* input, float* output, int n, cudaStream_t stream) {
				const size_t blockdim = next_pow_of_2(n / 2);	
				const size_t nblocks = (n + blockdim - 1) / blockdim;
				const size_t smem_size = blockdim * sizeof(float);
				const dim3 grid(nblocks, 1, 1);
				const dim3 blocks(nblocks, 1, 1);
        reduction_scalar<float><<<grid, blocks, smem_size, stream>>>(input, output, n); 
    }

    void reduction_scalar_double(double* input, double* output, int n, cudaStream_t stream) {
				const size_t blockdim = next_pow_of_2(n / 2);	
				const size_t nblocks = (n + blockdim - 1) / blockdim;
				const size_t smem_size = blockdim * sizeof(double);
				const dim3 grid(nblocks, 1, 1);
				const dim3 blocks(nblocks, 1, 1);
        reduction_scalar<double><<<grid, blocks, smem_size, stream>>>(input, output, n); 
    }

    void reduction_scalar_halo_float(float* input, float* output, int nh, int itot, int jtot, cudaStream_t stream) {
				const size_t blockdim = next_pow_of_2(itot / 2);
				const size_t smem_size = blockdim * sizeof(float);
        const dim3 grid(1, jtot, 1);
        const dim3 block(blockdim, 1, 1);
        reduction_scalar_halo<float><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot);
    }

    void reduction_scalar_halo_double(double* input, double* output, int nh, int itot, int jtot, cudaStream_t stream) {
				const size_t blockdim = next_pow_of_2(itot / 2);
				const size_t smem_size = blockdim * sizeof(double);
        const dim3 grid(1, jtot, 1);
        const dim3 block(blockdim, 1, 1);
				reduction_scalar_halo<double><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot);
    }

    void reduction_profile_halo_float(float* input, float* output, int nh, int itot, int jtot, int ktot, 
																			cudaStream_t stream) {
				const unsigned int blockdim = next_pow_of_2(itot / 2);
				const size_t smem_size = blockdim * sizeof(float);
        const dim3 grid(1, jtot, ktot);
        const dim3 block(blockdim, 1, 1);
        reduction_profile_halo<float><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
    }

    void reduction_profile_halo_double(double* input, double* output, int nh, int itot, int jtot, int ktot, 
																			 cudaStream_t stream) {
				const unsigned int blockdim = next_pow_of_2(itot / 2);
				const size_t smem_size = blockdim * sizeof(double);
        const dim3 grid(1, jtot, ktot);
        const dim3 block(blockdim, 1, 1);
				reduction_profile_halo<double><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
    }

}

#ifdef ENABLE_MAIN

int main() {

    const int itot = 1024;
    const int jtot = 1024;
		const int ktot = 80;
    const int nh = 1;

    const int size = (itot + 2*nh)*(jtot + 2*nh)*ktot;
		const int ijtot = (itot + 2 * nh) * (jtot + 2 * nh);
    const int bytes = size * sizeof(float);

		float* h_input = new float[size];
		float* h_output = new float[ktot];

		for (int k = 0; k < ktot; k++) {
    		for (int j = 0; j < jtot + 2 * nh; j++) {
    		    for (int i = 0; i < itot + 2 * nh; i++) {
    		        if (i >= nh && i < itot + nh && j >= nh && j < jtot + nh) {
    		            h_input[i + j*(itot+2*nh) + k*ijtot] = static_cast<float>(1);
    		        } else {
    		            h_input[i + j*(itot+2*nh) + k*ijtot] = static_cast<float>(-5000);
    		        }
    		    }
    		}
		}

		int sum { 0 };
    
		for (int k = 0; k < ktot; k++) {
    		for (int j = 0; j < jtot; j++) {
						for (int i = 0; i < itot; i++) {
	  		  			size_t ii = nh * (itot + 2 * nh) + nh + i + j * (itot + 2 * nh);
	  		  			sum += h_input[ii];
						}
    		}
				std::cout << "k = " << k << ", sum = " << sum << "\n";
				sum = 0;
		}

    float* dev_input;
    float* dev_output;

    cudaMalloc(&dev_input, bytes);
    cudaMalloc(&dev_output, ktot*sizeof(float));

    cudaMemcpy(dev_input, h_input, bytes, cudaMemcpyHostToDevice);

		const size_t nthreads = max(16, min(next_pow_of_2(itot/2), 1024));

		std::cout << nthreads << "\n";

		const dim3 grid(1, jtot, ktot);
		const dim3 block(nthreads, 1, 1);
		const size_t smem_size = nthreads * sizeof(float);

		reduction_profile_halo<float><<<grid, block, smem_size>>>(dev_input, dev_output, nh, itot, jtot, ktot);

		for (int k = 0; k < ktot; k++) {
				h_output[k] = static_cast<float>(0);
		}

    cudaMemcpy(h_output, dev_output, ktot*sizeof(float), cudaMemcpyDeviceToHost);

		std::cout << "==== Device result ====" << "\n";
		for (int k = 0; k < ktot; k++) {
				std::cout << "k = " << k << ", sum = " << h_output[k] << "\n";
		}

    delete[] h_input;
    delete[] h_output;
    cudaFree(dev_input);
    cudaFree(dev_output);

    return 0;
}

#endif
