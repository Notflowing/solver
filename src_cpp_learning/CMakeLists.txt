cmake_minimum_required(VERSION 3.18)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CUDA_ARCHITECTURES 80)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin)

project(cudaPossionSolver LANGUAGES CXX CUDA)

set(exec main)

add_executable(${exec} CudaPossion.cu)
target_include_directories(${exec} PUBLIC . ../inc)
target_compile_options(${exec} PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)
target_compile_options(${exec} PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)

find_package(CUDAToolkit 11.0 REQUIRED)
message(STATUS "Found CUDAToolkit ${CUDAToolkit_VERSION} at ${CUDAToolkit_LIBRARIES}")
target_link_libraries(${exec} PUBLIC CUDA::cublas)