#!/bin/bash

set -e

# 1: JACOBI
# 2: Gauss_Seidel
# 3: Sor

# ./bin/solver_CPU 1 101
./bin/solver_CPU 1 400 400
# ./bin/solver_CPU 1 101 101 101

# ./bin/solver_CPU 2 101
./bin/solver_CPU 2 400 400
# ./bin/solver_CPU 2 101 101 101

# ./bin/solver_CPU 3 101
./bin/solver_CPU 3 400 400
# ./bin/solver_CPU 3 101 101 101

# ./bin/solver_CPU 1 401 401 401