#include <iostream>
#include <cuda.h>

inline int next_pow_of_2(int x) { return pow(2, ceil(log(x) / log(2))); };

enum reduction_t {
		sum_t,
		max_t,
		min_t
};

template<typename T, reduction_t F>
__device__ T reduction(T x, T y) {
		T retval;
		switch(F) {
				case (sum_t):
						retval = x + y;
						break;
				case (max_t):
						retval = fmax(x, y);
						break;
				case (min_t):
						retval = fmin(x, y);
						break;
		}
		return retval;
}

template<typename T, unsigned int nthreads, reduction_t F>
__inline__ __device__ void block_reduce(volatile T* vals, const int tid) {

		if (nthreads >= 512 && tid < 256) { vals[tid] = reduction<T, F>(vals[tid], vals[tid + 256]); }

		__syncthreads();

		if (nthreads >= 256 && tid < 128) { vals[tid] = reduction<T, F>(vals[tid], vals[tid + 128]); }

		__syncthreads();

		if (nthreads >= 128 && tid < 64) { vals[tid] = reduction<T, F>(vals[tid], vals[tid + 64]); }

		__syncthreads();

		if (tid < 32) {
				T result = vals[tid];

				if (nthreads >= 64) { result = reduction<T, F>(result, vals[tid + 32]); }

				result = reduction<T, F>(result, __shfl_down_sync(0xffffffff, result, 16));
				result = reduction<T, F>(result, __shfl_down_sync(0xffffffff, result, 8));
				result = reduction<T, F>(result, __shfl_down_sync(0xffffffff, result, 4));
				result = reduction<T, F>(result, __shfl_down_sync(0xffffffff, result, 2));
				result = reduction<T, F>(result, __shfl_down_sync(0xffffffff, result, 1));

				if (tid == 0) { vals[tid] = result; }
		}
}

template<typename T, unsigned int nthreads>
__global__ void reduction_profile_halo(T* input, T* output, int nh, int itot, int jtot, int ktot) {
		int i = threadIdx.x;
		int j = blockIdx.y * 2;
		int k = blockIdx.z;
		int ijk = i + j * (itot + 2 * nh) + nh * (1 + itot + 2 * nh) + k * ((itot + 2 * nh) * (jtot + 2 * nh));
		
		extern __shared__ int smem_pointer[];
		T* smem = reinterpret_cast<T*>(smem_pointer);

		T sum = static_cast<T>(0);

		for (unsigned int tile = 0; tile < itot / blockDim.x + 1; tile++) {
				unsigned int offset = tile * blockDim.x;
				if (i + offset < itot) {
						sum += input[ijk + offset];
						if (j < jtot - 1) {
								sum += input[ijk + offset + itot + 2 * nh];
						}
				}	
		}

		smem[i] = sum;
		
		__syncthreads();

		block_reduce<T, nthreads, sum_t>(smem, i);

		if (i == 0) {
				atomicAdd(&output[k], smem[0]);
		}
}

// Interface part, which we link to the Fortran code
extern "C" {

    void reduction_profile_halo_float(float* input, float* output, int nh, int itot, int jtot, int ktot, 
																			cudaStream_t stream) {
				const unsigned int blockdim = next_pow_of_2(itot / 2);
				const size_t smem_size = blockdim * sizeof(float);
        const dim3 grid(1, jtot / 2, ktot);
        const dim3 block(blockdim, 1, 1);
				switch (blockdim) {
						case (512):
        				reduction_profile_halo<float, 512><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
						case (256):
        				reduction_profile_halo<float, 256><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
						case (128):
        				reduction_profile_halo<float, 128><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
						case (64):
        				reduction_profile_halo<float, 64><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
						case (32):
        				reduction_profile_halo<float, 32><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
						case (16):
        				reduction_profile_halo<float, 16><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
				}
    }

    void reduction_profile_halo_double(double* input, double* output, int nh, int itot, int jtot, int ktot, 
																			 cudaStream_t stream) {
				const unsigned int blockdim = next_pow_of_2(itot / 2);
				const size_t smem_size = blockdim * sizeof(double);
        const dim3 grid(1, jtot / 2, ktot);
        const dim3 block(blockdim, 1, 1);
				switch (blockdim) {
						case (512):
        				reduction_profile_halo<double, 512><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
						case (256):
        				reduction_profile_halo<double, 256><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
						case (128):
        				reduction_profile_halo<double, 128><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
						case (64):
        				reduction_profile_halo<double, 64><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
						case (32):
        				reduction_profile_halo<double, 32><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
						case (16):
        				reduction_profile_halo<double, 16><<<grid, block, smem_size, stream>>>(input, output, nh, itot, jtot, ktot);
								break;
				}
    }

}

int main() {

    constexpr int itot = 1024;
    constexpr int jtot = 1024;
		constexpr int ktot = 5;
    constexpr int nh = 1;

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

		const size_t nthreads = max(16, min(next_pow_of_2(itot/2), 512));

		std::cout << nthreads << "\n";

		const dim3 grid(1, jtot / 2, ktot);
		const dim3 block(nthreads, 1, 1);
		const size_t smem_size = nthreads * sizeof(float);

		reduction_profile_halo<float, itot / 2><<<grid, block, smem_size>>>(dev_input, dev_output, nh, itot, jtot, ktot);

		cudaDeviceSynchronize();

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
