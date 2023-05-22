#!/bin/bash

set -e


g++ -std=c++11 -DFP64 -O1 ./src_MG/main.cpp -I./inc -lm -o ./bin/main_MG
./bin/main_MG 5 400 400 | tee ./logDir/log_CPU_MG

# /public/software/cuda-11.5/bin/nvcc -std=c++11 -DFP32 -DGPU -arch=sm_80 -rdc=true -maxrregcount=127 -O2 ./src_MG/main_GPU.cpp -I./inc -I/public/software/cuda-11.5/include -lm -L/public/software/cuda-11.5/lib64 -lcudart -lcublas -o main_MG_GPU
# ./main_MG_GPU 5 400 400

make clean -f ./buildDir/Makefile_GPU_MG
make -f ./buildDir/Makefile_GPU_MG
bash ./buildDir/run_GPU_MG.sh | tee ./logDir/log_GPU_MG