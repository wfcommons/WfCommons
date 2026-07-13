#include "kernels.h"
#include "zero_one.h"
#include <hipcub/hipcub.hpp>

#include <chrono>
#include <iostream>
#include <optional>

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

// Function to run the GPU benchmark for a specified time
void runBenchmarkTime(long max_work, std::optional<int> runtime_in_seconds) {
  uint32_t n = 256 * 256;
  uint64_t m = (max_work + n - 1) / n;

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

  // set kernel
  dim3 gridSize = 256;
  dim3 blockSize = 256;

  setup_kernel<<<gridSize, blockSize>>>(d_state);

  HIP_CHECK(hipEventRecord(gpu_start, 0));
  int iteration = 0;

  if (!runtime_in_seconds.has_value()) {
    monte_carlo_kernel<<<gridSize, blockSize>>>(d_state, d_count, m);
    HIP_CHECK(hipDeviceSynchronize());
    iteration++;
  } else {
    // Run the workload loop until the specified runtime is reached
    while (getElapsedTime(gpu_start, gpu_stop) < runtime_in_seconds) {
      monte_carlo_kernel<<<gridSize, blockSize>>>(d_state, d_count, m);
      HIP_CHECK(hipDeviceSynchronize()); // Ensure the kernel has finished executing
      iteration++;
    }
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

void runBenchmark(long max_work) {
  runBenchmarkTime(max_work, std::nullopt);
}

// Function to run the GPU benchmark for a specified time
void memoryBenchmarkTime(long max_work, std::optional<int> runtime_in_seconds) {
  uint32_t n = 512;
  uint64_t m = (max_work + n - 1) / n;

  // allocate memory
  unsigned long long int *d_tmp;
  HIP_CHECK(hipMalloc((void **)&d_tmp, max_work * sizeof(unsigned long long int)));

  // set up timing stuff
  hipEvent_t gpu_start, gpu_stop;
  HIP_CHECK(hipEventCreate(&gpu_start));
  HIP_CHECK(hipEventCreate(&gpu_stop));

  // set kernel
  dim3 gridSize = m;
  dim3 blockSize = n;

  std::cout << n << " " << m << std::endl;

  // allocate memory
  unsigned long long int *d_count;
  HIP_CHECK(hipMalloc((void **)&d_count, m * sizeof(unsigned long long int)));
  HIP_CHECK(hipMemset(d_count, 0, m * sizeof(unsigned long long int)));

  // monte carlo kernel
  HIP_CHECK(hipEventRecord(gpu_start, 0));

  int iteration = 0;
  if (!runtime_in_seconds.has_value()) {
    HIP_CHECK(hipMemset(d_tmp, 0, max_work * sizeof(unsigned long long int)));
    stress_vm_zero<<<gridSize, blockSize>>>(d_tmp, d_count, max_work);
    HIP_CHECK(hipMemset(d_tmp, 0xff, max_work * sizeof(unsigned long long int)));
    stress_vm_one<<<gridSize, blockSize>>>(d_tmp, d_count, max_work);
    HIP_CHECK(hipDeviceSynchronize());
    iteration++;
  } else {
    // Run the workload loop until the specified runtime is reached
    while (getElapsedTime(gpu_start, gpu_stop) < runtime_in_seconds) {
      HIP_CHECK(hipMemset(d_tmp, 0, max_work * sizeof(unsigned long long int)));
      stress_vm_zero<<<gridSize, blockSize>>>(d_tmp, d_count, max_work);
      HIP_CHECK(hipMemset(d_tmp, 0xff, max_work * sizeof(unsigned long long int)));
      stress_vm_one<<<gridSize, blockSize>>>(d_tmp, d_count, max_work);
      HIP_CHECK(hipDeviceSynchronize()); // Ensure the kernel has finished executing
      iteration++;
    }
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
  HIP_CHECK(hipcub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_count, d_out, m));
  HIP_CHECK(hipMalloc((void **)&d_temp_storage, temp_storage_bytes));

  // Run
  HIP_CHECK(hipcub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, d_count, d_out, m));

  // copy results back to the host
  unsigned long long int h_count = 0;
  HIP_CHECK(hipMemcpy(&h_count, d_count, sizeof(unsigned long long int), hipMemcpyDeviceToHost));

  // display results and timings for gpu
  float pi = (float)h_count;
  std::cout << "Approximate pi calculated on GPU is: " << pi << " and calculation took " << gpu_elapsed_time << "s\n";

  HIP_CHECK(hipFree(d_tmp));
  HIP_CHECK(hipFree(d_count));
  HIP_CHECK(hipFree(d_temp_storage));
}

void memoryBenchmark(long max_work) {
  memoryBenchmarkTime(max_work, std::nullopt);
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

    memoryBenchmark(max_work);

  } else if (argc == 3) {
    // Parse the command line arguments
    long max_work = std::atol(argv[1]);
    int runtime_in_seconds = std::atoi(argv[2]);

    // Validate the input arguments
    if (max_work <= 0 || runtime_in_seconds <= 0) {
      std::cerr << "Both max_work and runtime_in_seconds must be positive integers." << std::endl;
      return 1;
    }

    memoryBenchmarkTime(max_work, runtime_in_seconds);

  } else {
    std::cerr << "Usage: " << argv[0] << " <max_work> [runtime_in_seconds]" << std::endl;
    return 1;
  }

  return 0;
}
