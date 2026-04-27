#include "kernels.h"

__global__ void setup_kernel(curandState *state) {
  int index = threadIdx.x + blockDim.x * blockIdx.x;
  curand_init(123456789, index, 0, &state[index]);
}

__global__ void monte_carlo_kernel(curandState *state, int64_t *count, int64_t m) {
  unsigned int index = threadIdx.x + blockDim.x * blockIdx.x;

  __shared__ int64_t cache[256];
  cache[threadIdx.x] = 0;
  __syncthreads();

  unsigned int temp = 0;
  while (temp < m) {
    float x = curand_uniform(&state[index]);
    float y = curand_uniform(&state[index]);
    float r = x * x + y * y;

    if (r <= 1) {
      cache[threadIdx.x]++;
    }
    temp++;
  }

  int i = blockDim.x / 2;
  while (i != 0) {
    if (threadIdx.x < i) {
      cache[threadIdx.x] += cache[threadIdx.x + i];
    }

    i /= 2;
    __syncthreads();
  }

  // update to our global variable count
  if (threadIdx.x == 0) {
    atomicAdd((unsigned long long int*)count, (unsigned long long int)cache[0]);
  }
}
