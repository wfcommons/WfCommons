#ifndef __ZERO_ONE_CUH__
#define __ZERO_ONE_CUH__

__global__ void stress_vm_zero(unsigned long long int *buf,
                               unsigned long long int *count,
                               const size_t sz);
__global__ void stress_vm_one(unsigned long long int *buf,
                              unsigned long long int *count,
                              const size_t sz);
#endif
