#ifndef __KERNELS_CUH__
#define __KERNELS_CUH__

#include <curand_kernel.h>

__global__ void setup_kernel(curandState *state);
__global__ void monte_carlo_kernel(curandState *state, unsigned long long int *count, int64_t m);

#endif
