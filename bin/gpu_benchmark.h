#ifndef GPU_BENCHMARK_H
#define GPU_BENCHMARK_H

#include <cuda_runtime.h>

void runBenchmark(int max_work);
void runBenchmarkTime(int max_work, int runtime_in_seconds);

#endif // GPU_BENCHMARK_H


