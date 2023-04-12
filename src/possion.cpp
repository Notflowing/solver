#include "header.h"

void init_u_host(type_t * __restrict__ u_cur, type_t * __restrict__ u_new,
                 const int Mdim, const int Ndim, const int Kdim)
{
    int i, j, k;
    long long idx;
    
    for (i = 0; i < Mdim; i++) {
        for (j = 0; j < Ndim; j++) {
            for (k = 0; k < Kdim; k++) {
                // idx = k + j * Kdim + i * Kdim * Ndim;
                u_cur[INDEX(i, j, k)] = 0.0;
                u_new[INDEX(i, j, k)] = 0.0;

            }
        }
    }
}

void initU(type_t *u_cur, type_t *u_new, const int Mdim, const int Ndim, const int Kdim)
{
// #if defined(__CUDA_ARCH__)
//     // init_u_kernel()
// #else
    init_u_host(u_cur, u_new, Mdim, Ndim, Kdim);
// #endif

}



void init_f_host(type_t * __restrict__ f, const int Mdim, const int Ndim, const int Kdim)
{
    int i, j, k;
    long long idx;
    
    for (i = 0; i < Mdim; i++) {
        for (j = 0; j < Ndim; j++) {
            for (k = 0; k < Kdim; k++) {
                idx = INDEX(i, j, k);
                f[idx] = 0.0;
                if (i == Mdim / 2 && j == Ndim / 2 && k == Kdim / 2) {
                    f[idx] = -(Kdim-1);
                }

            }
        }
    }
    
}

void initF(type_t *f, const int Mdim, const int Ndim, const int Kdim)
{
// #if defined(__CUDA_ARCH__)
//     // init_f_kernel()
// #else
    init_f_host(f, Mdim, Ndim, Kdim);
// #endif


}