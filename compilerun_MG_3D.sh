#!/bin/bash

set -e

# g++ -std=c++11 -DFP64 -O2 ./src_MG/main_3D.cpp -I./inc -lm -o ./bin/main_MG_3D
# ./bin/main_MG_3D 5 400 400 400 | tee ./logDir/log_CPU_MG_3D

make clean -f ./buildDir/Makefile_GPU_MG_3D
make -f ./buildDir/Makefile_GPU_MG_3D
bash ./buildDir/run_GPU_MG_3D.sh | tee ./logDir/log_GPU_MG_3D