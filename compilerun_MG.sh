#!/bin/bash

set -e


g++ -std=c++11 -DFP64 -O2 ./src_MG/main.cpp -I./inc -lm -o ./bin/main_MG
./bin/main_MG 5 1000 1000 | tee ./logDir/log_CPU_MG

make clean -f ./buildDir/Makefile_GPU_MG
make -f ./buildDir/Makefile_GPU_MG
bash ./buildDir/run_GPU_MG.sh | tee ./logDir/log_GPU_MG