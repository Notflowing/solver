#ifndef HEADER_H
#define HEADER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

// #include <string>
// #include <map>
// #include <tuple>
#include <vector>

#if defined(__CUDA_ARCH__)
    #include <cuda_runtime.h>
    #include "helper_cuda.h"
    #include "helper_string.h"
#else
    #define __host__  
    #define __device__  
    #define __forceinline__  
    #define __restrict__  
#endif


#if defined(FP64)
    #define type_t double
#elif defined(FP32)
    #define type_t float
#elif defined(FP16) && defined(__CUDA_ARCH__)
    #include "cuda_fp16.h"
    #define type_t half
#else
    #define type_t float
#endif

#define INDEX(i, j, k) ( (k) + (j) * (Kdim) + (i) * (Kdim) * (Ndim) )

#include "macro.h"
typedef void (*iterStencil_t)(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
                              const double DH2, const int Mdim, const int Ndim, const int Kdim,
                              const int i, const int j, const int k);
#include "function.h"
enum  ITERMET {JACOBI = 1, GAUSSSEIDEL = 2, SOR = 3, RBGS = 4, MGRBGS = 5};

#endif