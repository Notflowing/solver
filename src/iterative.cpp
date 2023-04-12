#include "header.h"

// Jacobi iteration method
__host__ __device__ __forceinline__
void jacobi_1D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
               const double DH2, const int Mdim, const int Ndim, const int Kdim,
               const int i, const int j, const int k)
{
    u_new[INDEX(i, j, k)] = 0.5 * ( u_cur[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] - DH2 * f[INDEX(i, j, k)] );
}

__host__ __device__ __forceinline__
void jacobi_2D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
               const double DH2, const int Mdim, const int Ndim, const int Kdim,
               const int i, const int j, const int k)
{
    u_new[INDEX(i, j, k)] = 0.25 * ( u_cur[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] + \
                                     u_cur[INDEX(i, j-1, k)] + u_cur[INDEX(i, j+1, k)] - DH2 * f[INDEX(i, j, k)] );
}

__host__ __device__ __forceinline__
void jacobi_3D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
               const double DH2, const int Mdim, const int Ndim, const int Kdim,
               const int i, const int j, const int k)
{
    u_new[INDEX(i, j, k)] = 1.0 / 6.0 * ( u_cur[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] + \
                                          u_cur[INDEX(i, j-1, k)] + u_cur[INDEX(i, j+1, k)] + \
                                          u_cur[INDEX(i-1, j, k)] + u_cur[INDEX(i+1, j, k)] - DH2 * f[INDEX(i, j, k)] );
}


// Gauss-Seidel iteration method
__host__ __device__ __forceinline__
void gaussseidel_1D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
                    const double DH2, const int Mdim, const int Ndim, const int Kdim,
                    const int i, const int j, const int k)
{
    u_new[INDEX(i, j, k)] = 0.5 * ( u_new[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] - DH2 * f[INDEX(i, j, k)] );
}

__host__ __device__ __forceinline__
void gaussseidel_2D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
                    const double DH2, const int Mdim, const int Ndim, const int Kdim,
                    const int i, const int j, const int k)
{
    u_new[INDEX(i, j, k)] = 0.25 * ( u_new[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] + \
                                     u_new[INDEX(i, j-1, k)] + u_cur[INDEX(i, j+1, k)] - DH2 * f[INDEX(i, j, k)] );
}

__host__ __device__ __forceinline__
void gaussseidel_3D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
               const double DH2, const int Mdim, const int Ndim, const int Kdim,
               const int i, const int j, const int k)
{
    u_new[INDEX(i, j, k)] = 1.0 / 6.0 * ( u_new[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] + \
                                          u_new[INDEX(i, j-1, k)] + u_cur[INDEX(i, j+1, k)] + \
                                          u_new[INDEX(i-1, j, k)] + u_cur[INDEX(i+1, j, k)] - DH2 * f[INDEX(i, j, k)] );
}


// SOR iteration method
__host__ __device__ __forceinline__
void sor_1D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
            const double DH2, const int Mdim, const int Ndim, const int Kdim,
            const int i, const int j, const int k)
{
    u_new[INDEX(i, j, k)] = RELAX * 0.5 * ( u_new[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] - \
                                            DH2 * f[INDEX(i, j, k)] ) + (1 - RELAX) * u_cur[INDEX(i, j, k)];
}

__host__ __device__ __forceinline__
void sor_2D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
            const double DH2, const int Mdim, const int Ndim, const int Kdim,
            const int i, const int j, const int k)
{
    u_new[INDEX(i, j, k)] = RELAX * 0.25 * ( u_new[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] + \
                                             u_new[INDEX(i, j-1, k)] + u_cur[INDEX(i, j+1, k)] - \
                                             DH2 * f[INDEX(i, j, k)] ) + (1 - RELAX) * u_cur[INDEX(i, j, k)];
}

__host__ __device__ __forceinline__
void sor_3D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
            const double DH2, const int Mdim, const int Ndim, const int Kdim,
            const int i, const int j, const int k)
{
    u_new[INDEX(i, j, k)] = RELAX * 1.0 / 6.0 * ( u_new[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] + \
                                                  u_new[INDEX(i, j-1, k)] + u_cur[INDEX(i, j+1, k)] + \
                                                  u_new[INDEX(i-1, j, k)] + u_cur[INDEX(i+1, j, k)] - \
                                                  DH2 * f[INDEX(i, j, k)] ) + (1 - RELAX) * u_cur[INDEX(i, j, k)];
}


void iteration_host(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
                    PARAMS params, INNERS inners, iterStencil_t iterstcl, double *delta)
{
    int Mdim = params.M;
    int Ndim = params.N;
    int Kdim = params.K;
    int DIM = params.DIM;
    double DH = params.DH;
    double DH2 = DH * DH;
    // iterStencil_t jacobistencil;
    *delta = 0.0;
    int istart = inners.istart;     int iend = inners.iend;
    int jstart = inners.jstart;     int jend = inners.jend;
    int kstart = inners.kstart;     int kend = inners.kend;

    int i, j, k;
    for (i = istart; i < iend; i++) {
        for (j = jstart; j < jend; j++) {
            for (k = kstart; k < kend; k++) {
                iterstcl(u_cur, u_new, f, DH2, Mdim, Ndim, Kdim, i, j, k);
                *delta += ( (u_new[INDEX(i, j, k)] - u_cur[INDEX(i, j, k)]) * (u_new[INDEX(i, j, k)] - u_cur[INDEX(i, j, k)]) );
            }
        }
    }
    *delta = sqrt((*delta));

}



double iterMethod(type_t *u_cur, type_t *u_new, type_t *f, PARAMS params,
                  INNERS inners, iterStencil_t iterstcl)
{
    double delta;

#if defined(__CUDA_ARCH__)
    iteration_dev (u_cur, u_new, f, params, inners, iterstcl, &delta);
#else
    iteration_host(u_cur, u_new, f, params, inners, iterstcl, &delta);
#endif

    return delta;

}
