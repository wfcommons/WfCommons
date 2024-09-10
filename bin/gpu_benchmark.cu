#include <iostream>
#include <chrono>
#include <cstdlib>  // For std::atoi
#include "gpu_benchmark.h"

// Kernel function to perform a simple workload
__global__ void simpleKernel(int* data, int size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        data[idx] = data[idx] * data[idx];  // Simple workload: squaring each element
    }
}

// Function to run the GPU benchmark
void runBenchmark(int max_work, int runtime_in_seconds) {
    int* h_data = new int[max_work];
    int* d_data;

    // Initialize data
    for (int i = 0; i < max_work; i++) {
        h_data[i] = i;
    }

    // Allocate GPU memory
    cudaMalloc(&d_data, max_work * sizeof(int));

    // Copy data to GPU
    cudaMemcpy(d_data, h_data, max_work * sizeof(int), cudaMemcpyHostToDevice);

    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    // Kernel configuration
    int threadsPerBlock = 256;
    int blocksPerGrid = (max_work + threadsPerBlock - 1) / threadsPerBlock;

    // Run the workload loop until the specified runtime is reached
    while (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() < runtime_in_seconds) {
        simpleKernel<<<blocksPerGrid, threadsPerBlock>>>(d_data, max_work);
        cudaDeviceSynchronize();  // Ensure the kernel has finished executing
    }

    // Copy results back to host (optional, just for validation)
    cudaMemcpy(h_data, d_data, max_work * sizeof(int), cudaMemcpyDeviceToHost);

    // Cleanup
    cudaFree(d_data);
    delete[] h_data;

    std::cout << "Benchmark completed!" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <max_work> <runtime_in_seconds>" << std::endl;
        return 1;
    }

    int max_work = std::atoi(argv[1]);
    int runtime_in_seconds = std::atoi(argv[2]);

    if (max_work <= 0 || runtime_in_seconds <= 0) {
        std::cerr << "Both max_work and runtime_in_seconds must be positive integers." << std::endl;
        return 1;
    }

    runBenchmark(max_work, runtime_in_seconds);

    return 0;
}
