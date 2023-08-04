#!/bin/bash

set -e

rm -r build
cmake -B build
cmake --build build
../bin/main
