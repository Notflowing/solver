#!/bin/bash

set -e

# 1: JACOBI
# 4: Red-Black Gauss_Seidel
# 5: MultiGrid with Red-Black Gauss_Seidel

# ./bin/solver_GPU 1 101
# ./bin/solver_GPU 1 101 101
# ./bin/solver_GPU 1 101 101 101

# ./bin/solver_GPU 4 101
# ./bin/solver_GPU 4 101 101
# ./bin/solver_GPU 4 101 101 101

# ./bin/solver_GPU 3 101
# ./bin/solver_GPU 3 101 101
# ./bin/solver_GPU 3 101 101 101

./bin/solver_GPU 1 401 401 401
./bin/solver_GPU 4 401 401 401