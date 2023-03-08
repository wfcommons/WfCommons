#include "kernels.cuh"




__global__ void setup_kernel(curandState *state)
{
	int index = threadIdx.x + blockDim.x*blockIdx.x;
    curand_init(123456789, index, 0, &state[index]);
}




__global__ void monte_carlo_kernel(curandState *state, int *count, int m)
{
	unsigned int index_x = threadIdx.x + blockDim.x*blockIdx.x;
	unsigned int index_y = threadIdx.y + blockDim.y*blockIdx.y;
	
	__shared__ int cache[16*16];
	// cache[threadIdx.x] = 0;
	// cache[threadIdx.y] = 0;

	// __syncthreads();

	
	unsigned int temp = 0;
	while(temp < m){
		unsigned u = threadIdx.y*blockDim.x + threadIdx.x;
		unsigned i = curand_uniform(&state[index_x]);
		unsigned j = curand_uniform(&state[index_y]);
		unsigned Ni = gridDim.x*blockDim.x;
		unsigned Nj = gridDim.y*blockDim.y;

		float x = i/(float)Ni;
		float y = j/(float)Nj;
		float r = std::sqrt(x*x + y*y);

		cache[u] += r<=1; 
		temp++; 
	}

	
	// reduction --- probably remove this part 
	// int i = blockDim.x/2;
	// while(i != 0){
	// 	if(threadIdx.x < i){
	// 		cache[threadIdx.x] += cache[threadIdx.x + i];
	// 	}

	// 	i /= 2;
	// 	__syncthreads();
	// }


	// // update to our global variable count
	// if(threadIdx.x == 0){
	// 	atomicAdd(count, cache[0]);
	// }
}



