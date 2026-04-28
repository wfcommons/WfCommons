#include "kernels.h"

#include <cub/cub.cuh>

__global__ void setup_kernel(curandState *state) {
  int index = threadIdx.x + blockDim.x * blockIdx.x;
  curand_init(123456789, index, 0, &state[index]);
}

__global__ void monte_carlo_kernel(curandState *state,
                                   unsigned long long int *count, int64_t m) {
  unsigned int index = threadIdx.x + blockDim.x * blockIdx.x;

  unsigned long long int thread_data = 0;

  unsigned int temp = 0;
  while (temp < m) {
    float x = curand_uniform(&state[index]);
    float y = curand_uniform(&state[index]);
    float r = x * x + y * y;

    if (r <= 1) {
      thread_data++;
    }
    temp++;
  }

  typedef cub::BlockReduce<unsigned long long int, 256> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp_storage;
  unsigned long long int aggregate =
      BlockReduceT(temp_storage).Sum(thread_data);

  // update to our global variable count
  if (threadIdx.x == 0) {
    count[blockIdx.x] += aggregate;
  }
}