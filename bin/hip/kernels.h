#include "hip/hip_runtime.h"
#ifndef __KERNELS_CUH__
#define __KERNELS_CUH__

#include <hiprand/hiprand_kernel.h>

__global__ void setup_kernel(hiprandState *state);
__global__ void monte_carlo_kernel(hiprandState *state, unsigned long long int *count, int64_t m);

#endif
