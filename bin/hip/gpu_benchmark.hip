#include "kernels.h"
#include <hipcub/hipcub.hpp>

#include <chrono>
#include <iostream>

#define HIP_CHECK(expression)                  \
{                                              \
    const hipError_t status = expression;      \
    if(status != hipSuccess){                  \
        std::cerr << "HIP error "              \
                  << status << ": "            \
                  << hipGetErrorString(status) \
                  << " at " << __FILE__ << ":" \
                  << __LINE__ << std::endl;    \
    }                                          \
}

float getElapsedTime(const hipEvent_t &gpu_start, hipEvent_t &gpu_stop) {
  float gpu_elapsed_time;
  HIP_CHECK(hipEventRecord(gpu_stop, 0));
  HIP_CHECK(hipEventSynchronize(gpu_stop));
  HIP_CHECK(hipEventElapsedTime(&gpu_elapsed_time, gpu_start, gpu_stop));
  return gpu_elapsed_time / 1000.0f;
}

// Function to run the GPU benchmark with no time limit
void runBenchmark(long max_work) {
  uint32_t n = 256 * 256;
  uint64_t m = max_work * 16384 / n;

  unsigned long long int *d_count;
  hiprandState *d_state;
  HIP_CHECK(hipMalloc((void **)&d_count, 256 * sizeof(unsigned long long int)));
  HIP_CHECK(hipMalloc((void **)&d_state, n * sizeof(hiprandState)));
  HIP_CHECK(hipMemset(d_count, 0, 256 * sizeof(unsigned long long int)));

  // set up timing stuff
  hipEvent_t gpu_start, gpu_stop;
  HIP_CHECK(hipEventCreate(&gpu_start));
  HIP_CHECK(hipEventCreate(&gpu_stop));

  // set kernel
  dim3 gridSize = 256;
  dim3 blockSize = 256;
  setup_kernel<<<gridSize, blockSize>>>(d_state);

  // monte carlo kernel
  HIP_CHECK(hipEventRecord(gpu_start, 0));
  monte_carlo_kernel<<<gridSize, blockSize>>>(d_state, d_count, m);
  HIP_CHECK(hipDeviceSynchronize());

  float gpu_elapsed_time = getElapsedTime(gpu_start, gpu_stop);
  HIP_CHECK(hipEventDestroy(gpu_start));
  HIP_CHECK(hipEventDestroy(gpu_stop));

  // Allocate device output array
  unsigned long long int *d_out = nullptr;
  HIP_CHECK(hipMalloc((void **)&d_out, sizeof(unsigned long long int)));

  // Request and allocate temporary storage
  void *d_temp_storage = nullptr;
  size_t temp_storage_bytes = 0;
  HIP_CHECK(hipcub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_count, d_out, 256));
  HIP_CHECK(hipMalloc((void **)&d_temp_storage, temp_storage_bytes));

  // Run
  HIP_CHECK(hipcub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_count, d_out, 256));

  // copy results back to the host
  unsigned long long int h_count = 0;
  HIP_CHECK(hipMemcpy(&h_count, d_out, sizeof(unsigned long long int), hipMemcpyDeviceToHost));

  // display results and timings for gpu
  float pi = h_count * 4.0 / (n * m);
  std::cout << "Approximate pi calculated on GPU is: " << pi << " and calculation took " << gpu_elapsed_time << "s\n";
  std::cout << "Benchmark completed!" << std::endl;

  HIP_CHECK(hipFree(d_count));
  HIP_CHECK(hipFree(d_state));
  HIP_CHECK(hipFree(d_out));
  HIP_CHECK(hipFree(d_temp_storage));
}

// Function to run the GPU benchmark for a specified time
void runBenchmarkTime(long max_work, int runtime_in_seconds) {

  uint32_t n = 256 * 256;
  uint64_t m = max_work * 16384 / n;

  // allocate memory
  unsigned long long int *d_count;
  hiprandState *d_state;
  HIP_CHECK(hipMalloc((void **)&d_count, 256 * sizeof(unsigned long long int)));
  HIP_CHECK(hipMalloc((void **)&d_state, n * sizeof(hiprandState)));
  HIP_CHECK(hipMemset(d_count, 0, 256 * sizeof(unsigned long long int)));

  // set up timing stuff
  hipEvent_t gpu_start, gpu_stop;
  HIP_CHECK(hipEventCreate(&gpu_start));
  HIP_CHECK(hipEventCreate(&gpu_stop));
  HIP_CHECK(hipEventRecord(gpu_start, 0));

  // set kernel
  dim3 gridSize = 256;
  dim3 blockSize = 256;

  setup_kernel<<<gridSize, blockSize>>>(d_state);

  int iteration = 0;
  // Run the workload loop until the specified runtime is reached
  while (getElapsedTime(gpu_start, gpu_stop) < runtime_in_seconds) {
    monte_carlo_kernel<<<gridSize, blockSize>>>(d_state, d_count, m);
    HIP_CHECK(hipDeviceSynchronize()); // Ensure the kernel has finished executing
    iteration++;
  }

  float gpu_elapsed_time = getElapsedTime(gpu_start, gpu_stop);
  HIP_CHECK(hipEventDestroy(gpu_start));
  HIP_CHECK(hipEventDestroy(gpu_stop));

  // copy results back to the host
  // Allocate device output array
  unsigned long long int *d_out = nullptr;
  HIP_CHECK(hipMalloc((void **)&d_out, sizeof(unsigned long long int)));

  // Request and allocate temporary storage
  void *d_temp_storage = nullptr;
  size_t temp_storage_bytes = 0;
  HIP_CHECK(hipcub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_count, d_out, 256));
  HIP_CHECK(hipMalloc((void **)&d_temp_storage, temp_storage_bytes));

  // Run
  HIP_CHECK(hipcub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_count, d_out, 256));

  // copy results back to the host
  unsigned long long int h_count = 0;
  HIP_CHECK(hipMemcpy(&h_count, d_out, sizeof(unsigned long long int), hipMemcpyDeviceToHost));

  // display results and timings for gpu
  float pi = h_count * 4.0 / (n * m) / iteration;
  std::cout << "Approximate pi calculated on GPU is: " << pi << " and calculation took " << gpu_elapsed_time << "s\n";

  HIP_CHECK(hipFree(d_count));
  HIP_CHECK(hipFree(d_state));
  HIP_CHECK(hipFree(d_out));
  HIP_CHECK(hipFree(d_temp_storage));
}

int main(int argc, char *argv[]) {
  // Check for the correct number of command line arguments
  if (argc == 2) {
    // Parse the command line arguments
    long max_work = std::atol(argv[1]);

    // Validate the input arguments
    if (max_work <= 0) {
      std::cerr << "max_work must be a positive integer." << std::endl;
      return 1;
    }

    runBenchmark(max_work);

  } else if (argc == 3) {
    // Parse the command line arguments
    long max_work = std::atol(argv[1]);
    int runtime_in_seconds = std::atoi(argv[2]);

    // Validate the input arguments
    if (max_work <= 0 || runtime_in_seconds <= 0) {
      std::cerr << "Both max_work and runtime_in_seconds must be positive integers." << std::endl;
      return 1;
    }

    runBenchmarkTime(max_work, runtime_in_seconds);

  } else {
    std::cerr << "Usage: " << argv[0] << " <max_work> [runtime_in_seconds]" << std::endl;
    return 1;
  }

  return 0;
}
