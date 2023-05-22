#!/bin/bash

set -e

make clean -f ./buildDir/Makefile_GPU
make -f ./buildDir/Makefile_GPU
bash ./buildDir/run_GPU.sh | tee ./logDir/log_GPU