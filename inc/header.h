#ifndef HEADER_H
#define HEADER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

// #include <string>
// #include <map>
// #include <tuple>
#include <vector>
#include <functional>

#if defined(GPU)
    #include <cuda_runtime.h>
    #include <cublas_v2.h>
    // #include <thrust/universal_vector.h>
    #include "helper_cuda.h"
    #include "helper_string.h"
#endif

#if defined(FP64)
    #define type_t double
#elif defined(FP32)
    #define type_t float
#elif defined(FP16) && defined(GPU)
    #include "cuda_fp16.h"
    #define type_t half
#else
    #define type_t float
#endif

#define INDEX(i, j, k) ( (k) + (j) * (Kdim) + (i) * (Kdim) * (Ndim) )
#define INDEX_fine(i, j, k) ( (k) + (j) * (Kdim_fine) + (i) * (Kdim_fine) * (Ndim_fine) )
#define INDEX_coar(i, j, k) ( (k) + (j) * (Kdim_coar) + (i) * (Kdim_coar) * (Ndim_coar) )

#include "macro.h"
typedef void (*iterStencil_t)(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
                              const double DH2, const int Mdim, const int Ndim, const int Kdim,
                              const int i, const int j, const int k);
enum  ITERMET {JACOBI = 1, GAUSSSEIDEL = 2, SOR = 3, RBGS = 4, MGRBGS = 5};

#if defined(GPU)
extern __constant__ int dim_dev[3];
extern __constant__ int inn_dev[6];
extern cublasHandle_t handle;
#endif

#include "function.h"

#define KHZ 4
#define PI 3.1415926535


#endif