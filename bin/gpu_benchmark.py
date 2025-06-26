import numpy as np
import time
import warnings
import argparse
from numba import cuda, float32
from numba.core.errors import NumbaPerformanceWarning

warnings.simplefilter('ignore', category=NumbaPerformanceWarning)

def setup_device():
    """
    Prints the name of the currently active CUDA GPU device.
    """
    print("Using GPU:", cuda.get_current_device().name)


# ---------------------------------------------
# 1. Compute-bound: Matrix Multiplication
# ---------------------------------------------
@cuda.jit
def matmul_kernel(A, B, C):
    """
    CUDA kernel to perform matrix multiplication C = A * B.

    Each thread computes one element of the output matrix C.
    
    Args:
        A (device array): Left matrix operand with shape (N, N).
        B (device array): Right matrix operand with shape (N, N).
        C (device array): Output matrix to store the product with shape (N, N).
    """
    row, col = cuda.grid(2)
    if row < C.shape[0] and col < C.shape[1]:
        tmp = 0.
        for k in range(A.shape[1]):
            tmp += A[row, k] * B[k, col]
        C[row, col] = tmp

def benchmark_compute_bound(N=1024):
    """
    Benchmark measuring compute-bound GPU performance using matrix multiplication.

    Initializes two NxN random matrices, copies them to the GPU, runs
    matrix multiplication kernel twice (once for warmup, once timed),
    and prints the elapsed GPU execution time.

    Args:
        N (int): Size of the square matrices to multiply. Default is 1024.
    """
    A = np.random.rand(N, N).astype(np.float32)
    B = np.random.rand(N, N).astype(np.float32)

    d_A = cuda.to_device(A)
    d_B = cuda.to_device(B)
    d_C = cuda.device_array((N, N), dtype=np.float32)

    threads_per_block = (16, 16)
    blocks_per_grid = (N // 16, N // 16)

    # Warm-up kernel launch to mitigate startup overhead
    matmul_kernel[blocks_per_grid, threads_per_block](d_A, d_B, d_C)
    cuda.synchronize()

    # Timed kernel launch
    start = time.time()
    matmul_kernel[blocks_per_grid, threads_per_block](d_A, d_B, d_C)
    cuda.synchronize()
    end = time.time()

    print(f"[Compute] Matrix multiplication time: {end - start:.4f} sec")


# ---------------------------------------------
# 2. Memory-bound: Bandwidth Test
# ---------------------------------------------
@cuda.jit
def memory_copy_kernel(src, dst):
    """
    CUDA kernel to copy elements from src array to dst array.

    Args:
        src (device array): Source array of floats.
        dst (device array): Destination array to copy data into.
    """
    i = cuda.grid(1)
    if i < src.size:
        dst[i] = src[i]

def benchmark_memory_bound(N=100_000_000):
    """
    Benchmark measuring memory bandwidth by copying large arrays on the GPU.

    Allocates two large arrays on host and device, runs a kernel to copy data
    from src to dst, and computes bandwidth in GB/s.

    Args:
        N (int): Number of float elements in the arrays. Default is 100 million.
    """
    src = np.random.rand(N).astype(np.float32)
    dst = np.zeros_like(src)

    d_src = cuda.to_device(src)
    d_dst = cuda.device_array_like(dst)

    threads_per_block = 256
    blocks_per_grid = (N + threads_per_block - 1) // threads_per_block

    # Warm-up kernel launch
    memory_copy_kernel[blocks_per_grid, threads_per_block](d_src, d_dst)
    cuda.synchronize()

    # Timed kernel launch
    start = time.time()
    memory_copy_kernel[blocks_per_grid, threads_per_block](d_src, d_dst)
    cuda.synchronize()
    end = time.time()

    bytes_moved = src.nbytes
    bandwidth = bytes_moved / (end - start) / 1e9  # Convert bytes/sec to GB/s
    print(f"[GPU Memory] Memory bandwidth: {bandwidth:.2f} GB/s")


# ---------------------------------------------
# 3. Kernel Launch Overhead
# ---------------------------------------------
@cuda.jit
def noop_kernel():
    """
    A minimal CUDA kernel that does no computation.
    Used to benchmark kernel launch overhead.
    """
    pass

def benchmark_kernel_overhead(repeats=10000, blocks=1, threads=1):
    """
    Benchmark measuring the average overhead of launching a minimal CUDA kernel.

    Launches the noop_kernel repeatedly and measures total and average launch times.

    Args:
        repeats (int): Number of kernel launches to time. Default is 10,000.
    """
    start = time.time()
    for _ in range(repeats):
        noop_kernel[blocks, threads]()  # Single thread block with one thread
    cuda.synchronize()
    end = time.time()

    avg_launch_time = (end - start) / repeats * 1e6  # microseconds per launch
    print(f"[GPU Launch] Average kernel launch overhead: {avg_launch_time:.2f} µs")


# ---------------------------------------------
# 4. Workflow-style Composite Task
# ---------------------------------------------
@cuda.jit
def scale_and_square_kernel(data):
    """
    CUDA kernel to scale and square each element of the input array in-place.

    Args:
        data (device array): Array of float32 elements.
    """
    i = cuda.grid(1)
    if i < data.size:
        val = data[i]
        data[i] = (val * 2.0) ** 2


# ---------------------------------------------
# Run all benchmarks
# ---------------------------------------------
def run_all_benchmarks():
    """
    Runs all GPU benchmark tests sequentially:
    1. Compute-bound matrix multiplication
    2. Memory bandwidth test
    3. Kernel launch overhead
    """
    benchmark_compute_bound()
    benchmark_memory_bound()
    benchmark_kernel_overhead()


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gpu_compute_size", help="Activate gpu compute benchmark")
    parser.add_argument("--gpu_mem_size", help="Activate gpu memory bandwidth benchmark")
    parser.add_argument("--gpu_kernel", action="store_true", help="Activate gpu kernel launch benchmark")
    parser.add_argument("--gpu_all", action="store_true", help="Activate all gpu benchmarks")
    parser.add_argument(
        "--gpu_kernel_size",
        nargs=2,
        type=int,
        metavar=('BLOCKS', 'THREADS'),
        default=[1, 1],
        help="Grid and block size for kernel overhead benchmark (default: 1 1)"
    )
    return parser

if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()
    comp_test = int(args.gpu_compute_size)
    mem_test = args.gpu_mem
    kernel_test = args.gpu_kernel
    kernel_size = args.gpu_kernel_size
    all_test = args.gpu_all

    
    setup_device()
    if comp_test:
        benchmark_compute_bound(comp_test)
    elif mem_test:
        benchmark_memory_bound(mem_test)
    elif kernel_test or kernel_size:
        blocks, threads = kernel_size
        benchmark_kernel_overhead(blocks=blocks, threads=threads)
    elif all_test:
        run_all_benchmarks()
