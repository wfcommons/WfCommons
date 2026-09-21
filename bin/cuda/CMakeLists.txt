cmake_minimum_required(VERSION 3.21) # HIP language support requires 3.21
cmake_policy(VERSION 3.21.3...3.27)

find_package(CUDAToolkit REQUIRED)

project(MyProj LANGUAGES CUDA)
add_executable(gpu_benchmark gpu_benchmark.cu kernels.cu zero_one.cu)
target_link_libraries(gpu_benchmark PRIVATE CUDA::cudart CUDA::cuda_driver) 
