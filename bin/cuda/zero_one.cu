#include "zero_one.h"

#include <cub/cub.cuh>

using uint_type = unsigned long long int;

/*
 *  stress_vm_zero()
 *	set all memory to zero and see if any bits are stuck at one and
 */
__global__ void stress_vm_zero(
	uint_type *buf,
	uint_type *count, 
	const size_t sz)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= sz)
        return;

	uint_type bit_errors = __popcll(buf[tid]);

	typedef cub::BlockReduce<unsigned long long int, 512> BlockReduceT;
    __shared__ typename BlockReduceT::TempStorage temp_storage;
    uint_type aggregate = BlockReduceT(temp_storage).Sum(bit_errors);

    // update to our global variable count
    if (threadIdx.x == 0) {
      count[blockIdx.x] += aggregate;
	}
}

/*
 *  stress_vm_one()
 *	set all memory to one and see if any bits are stuck at zero
 */
__global__ void stress_vm_one(
	uint_type *buf,
	uint_type *count, 
	const size_t sz)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid >= sz)
        return;

	uint_type bit_errors = __popcll(~buf[tid]);

	typedef cub::BlockReduce<unsigned long long int, 512> BlockReduceT;
    __shared__ typename BlockReduceT::TempStorage temp_storage;
    uint_type aggregate = BlockReduceT(temp_storage).Sum(bit_errors);

    // update to our global variable count
    if (threadIdx.x == 0) {
      count[blockIdx.x] += aggregate;
	}
}
