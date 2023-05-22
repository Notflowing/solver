#!/bin/bash

set -e

make clean -f ./buildDir/Makefile_CPU
make -f ./buildDir/Makefile_CPU
bash ./buildDir/run_CPU.sh | tee ./logDir/log_CPU