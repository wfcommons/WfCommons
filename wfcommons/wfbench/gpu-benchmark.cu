#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <curand.h>
#include <time.h>
#include <math.h>
#include "kernels.cuh"


int main(int argc, char** argv)
{	

	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <work (# 1M samples)>\n";
		exit(1);
    }

	unsigned int n = 256*256;
	unsigned int m;
	unsigned int work;
	int *h_count;
	int *d_count;
	curandState *d_state;
	// float pi;

	//checking the user input for work
	try {
    	work = std::stol(argv[1]);
    } catch (std::invalid_argument &e) {
        std::cerr << "Invalid argument: " << e.what() << "\n";
        exit(1);
    }

	//making into M samples
	m = 1000000*work;
	// allocate memory
	h_count = (int*)malloc(n*sizeof(int));
	cudaMalloc((void**)&d_count, n*sizeof(int));
	cudaMalloc((void**)&d_state, n*sizeof(curandState));
	cudaMemset(d_count, 0, sizeof(int));


	// set up timing stuff
	float gpu_elapsed_time;
	cudaEvent_t gpu_start, gpu_stop;
	cudaEventCreate(&gpu_start);
	cudaEventCreate(&gpu_stop);
	cudaEventRecord(gpu_start, 0);


	// set kernel
	dim3 gridSize = 16;
	dim3 blockSize = 16;
	setup_kernel<<< gridSize, blockSize>>>(d_state);


	// monti carlo kernel
	monte_carlo_kernel<<<gridSize, blockSize>>>(d_state, d_count, m);


	// // copy results back to the host
	// cudaMemcpy(h_count, d_count, sizeof(int), cudaMemcpyDeviceToHost);
	cudaEventRecord(gpu_stop, 0);
	cudaEventSynchronize(gpu_stop);
	cudaEventElapsedTime(&gpu_elapsed_time, gpu_start, gpu_stop);
	cudaEventDestroy(gpu_start);
	cudaEventDestroy(gpu_stop);


	// // display results and timings for gpu
	// pi = *h_count*4.0/(n*m);
	// std::cout<<"Approximate pi calculated on GPU is: "<<pi<<" and calculation took "<<gpu_elapsed_time<<std::endl;
	std::cout<<"GPU stress test is over and it took "<<gpu_elapsed_time<<std::endl;

	// delete memory
	free(h_count);
	cudaFree(d_count);
	cudaFree(d_state);
}

