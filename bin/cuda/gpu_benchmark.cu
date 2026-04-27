#include "gpu_benchmark.h"
#include "kernels.h"
#include <chrono>
#include <cstdlib> // For std::atoi
#include <iostream>

// The macro wraps any CUDA API call
#define CUDA_CHECK(ans)                                                        \
  { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line,
                      bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
            line);
    if (abort)
      exit(code);
  }
}

// Function to run the GPU benchmark with no time limit
void runBenchmark(int max_work) {
  uint32_t n = 256 * 256;
  uint64_t m = max_work * 16384 / n;
  int64_t *h_count;
  int64_t *d_count;
  curandState *d_state;
  float pi;

  // allocate memory
  h_count = (int64_t *)malloc(n * sizeof(int64_t));
  cudaMalloc((void **)&d_count, n * sizeof(int64_t));
  cudaMalloc((void **)&d_state, n * sizeof(curandState));
  cudaMemset(d_count, 0, sizeof(int64_t));

  // set up timing stuff
  float gpu_elapsed_time;
  cudaEvent_t gpu_start, gpu_stop;
  cudaEventCreate(&gpu_start);
  cudaEventCreate(&gpu_stop);
  cudaEventRecord(gpu_start, 0);

  // set kernel
  dim3 gridSize = 256;
  dim3 blockSize = 256;
  setup_kernel<<<gridSize, blockSize>>>(d_state);

  // monte carlo kernel
  monte_carlo_kernel<<<gridSize, blockSize>>>(d_state, d_count, m);
  cudaDeviceSynchronize();

  // copy results back to the host
  cudaMemcpy(h_count, d_count, sizeof(int64_t), cudaMemcpyDeviceToHost);
  cudaEventRecord(gpu_stop, 0);
  cudaEventSynchronize(gpu_stop);
  cudaEventElapsedTime(&gpu_elapsed_time, gpu_start, gpu_stop);
  cudaEventDestroy(gpu_start);
  cudaEventDestroy(gpu_stop);

  // display results and timings for gpu
  pi = *h_count * 4.0 / (n * m);
  std::cout << "Approximate pi calculated on GPU is: " << pi << " " << *h_count
            << " and calculation took " << gpu_elapsed_time << std::endl;

  std::cout << "Benchmark completed!" << std::endl;
}

// Function to run the GPU benchmark for a specified time
void runBenchmarkTime(int max_work, int runtime_in_seconds) {

  uint32_t n = 256 * 256;
  uint64_t m = max_work * 16384 / n;
  int64_t *h_count;
  int64_t *d_count;
  curandState *d_state;
  float pi;

  // allocate memory
  h_count = (int64_t *)malloc(n * sizeof(int64_t));
  cudaMalloc((void **)&d_count, n * sizeof(int64_t));
  cudaMalloc((void **)&d_state, n * sizeof(curandState));
  cudaMemset(d_count, 0, sizeof(int64_t));

  // set up timing stuff
  float gpu_elapsed_time;
  cudaEvent_t gpu_start, gpu_stop;
  cudaEventCreate(&gpu_start);
  cudaEventCreate(&gpu_stop);
  cudaEventRecord(gpu_start, 0);

  // set kernel
  dim3 gridSize = 256;
  dim3 blockSize = 256;

  auto start = std::chrono::high_resolution_clock::now();
  setup_kernel<<<gridSize, blockSize>>>(d_state);

  int iteration = 0;
  // Run the workload loop until the specified runtime is reached
  while (std::chrono::duration_cast<std::chrono::seconds>(
             std::chrono::high_resolution_clock::now() - start)
             .count() < runtime_in_seconds) {
    monte_carlo_kernel<<<gridSize, blockSize>>>(d_state, d_count, m);
    cudaDeviceSynchronize(); // Ensure the kernel has finished executing
    iteration++;
  }

  // copy results back to the host
  cudaMemcpy(h_count, d_count, sizeof(int64_t), cudaMemcpyDeviceToHost);
  cudaEventRecord(gpu_stop, 0);
  cudaEventSynchronize(gpu_stop);
  cudaEventElapsedTime(&gpu_elapsed_time, gpu_start, gpu_stop);
  cudaEventDestroy(gpu_start);
  cudaEventDestroy(gpu_stop);

  // display results and timings for gpu
  pi = *h_count * 4.0 / (n * m) / iteration;
  std::cout << "Approximate pi calculated on GPU is: " << pi
            << " and calculation took " << gpu_elapsed_time << std::endl;

  std::cout << "Benchmark completed!" << std::endl;
}

int main(int argc, char *argv[]) {
  // Check for the correct number of command line arguments
  if (argc == 2) {
    // Parse the command line arguments
    int max_work = std::atoi(argv[1]);

    // Validate the input arguments
    if (max_work <= 0) {
      std::cerr << "max_work must be a positive integer." << std::endl;
      return 1;
    }

    runBenchmark(max_work);

  } else if (argc == 3) {
    // Parse the command line arguments
    int max_work = std::atoi(argv[1]);
    int runtime_in_seconds = std::atoi(argv[2]);

    // Validate the input arguments
    if (max_work <= 0 || runtime_in_seconds <= 0) {
      std::cerr
          << "Both max_work and runtime_in_seconds must be positive integers."
          << std::endl;
      return 1;
    }

    runBenchmarkTime(max_work, runtime_in_seconds);

  } else {
    std::cerr << "Usage: " << argv[0] << " <max_work> [runtime_in_seconds]"
              << std::endl;
    return 1;
  }

  return 0;
}
