#include "header.h"


__global__ void jacobi_1D_kernel(type_t * __restrict__ u_cur, type_t * __restrict__ u_new,
                                 type_t * __restrict__ error, type_t * __restrict__ f,
                                 const double DH2)
{
    int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int i = threadIdx.z + blockIdx.z * blockDim.z;

    int Mdim = dim_dev[0];
    int Ndim = dim_dev[1];
    int Kdim = dim_dev[2];

    if (i >= inn_dev[0] && i < inn_dev[1] && j >= inn_dev[2] && j < inn_dev[3] && k >= inn_dev[4] && k < inn_dev[5]) {
        // printf("i=%d, j=%d, k=%d\n", i, j, k);
        u_new[INDEX(i, j, k)] = 0.5 * ( u_cur[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] - DH2 * f[INDEX(i, j, k)] );
        error[INDEX(i, j, k)] = u_new[INDEX(i, j, k)] - u_cur[INDEX(i, j, k)];
    }

}


__global__ void jacobi_2D_kernel(type_t * __restrict__ u_cur, type_t * __restrict__ u_new,
                                 type_t * __restrict__ error, type_t * __restrict__ f,
                                 const double DH2)
{
    int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
    int j = threadIdx.y + blockIdx.y * blockDim.y + 1;
    int i = threadIdx.z + blockIdx.z * blockDim.z;

    int Mdim = dim_dev[0];
    int Ndim = dim_dev[1];
    int Kdim = dim_dev[2];

    if (i >= inn_dev[0] && i < inn_dev[1] && j >= inn_dev[2] && j < inn_dev[3] && k >= inn_dev[4] && k < inn_dev[5]) {
        u_new[INDEX(i, j, k)] = 0.25 * ( u_cur[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] + \
                                         u_cur[INDEX(i, j-1, k)] + u_cur[INDEX(i, j+1, k)] - DH2 * f[INDEX(i, j, k)] );
        error[INDEX(i, j, k)] = u_new[INDEX(i, j, k)] - u_cur[INDEX(i, j, k)];
    }

}

__global__ void jacobi_3D_kernel(type_t * __restrict__ u_cur, type_t * __restrict__ u_new,
                                 type_t * __restrict__ error, type_t * __restrict__ f,
                                 const double DH2)
{
    int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
    int j = threadIdx.y + blockIdx.y * blockDim.y + 1;
    int i = threadIdx.z + blockIdx.z * blockDim.z + 1;

    int Mdim = dim_dev[0];
    int Ndim = dim_dev[1];
    int Kdim = dim_dev[2];

    if (i >= inn_dev[0] && i < inn_dev[1] && j >= inn_dev[2] && j < inn_dev[3] && k >= inn_dev[4] && k < inn_dev[5]) {
        u_new[INDEX(i, j, k)] = 1.0 / 6.0 * ( u_cur[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] + \
                                              u_cur[INDEX(i, j-1, k)] + u_cur[INDEX(i, j+1, k)] + \
                                              u_cur[INDEX(i-1, j, k)] + u_cur[INDEX(i+1, j, k)] - DH2 * f[INDEX(i, j, k)] );
        error[INDEX(i, j, k)] = u_new[INDEX(i, j, k)] - u_cur[INDEX(i, j, k)];
    }

}


template <int phase>
__global__ void rbgs_1D_kernel(type_t * __restrict__ u_cur, type_t * __restrict__ f, const double DH2)
{
    // int k = (threadIdx.x + blockIdx.x * blockDim.x) * 2 + phase + 1;
    int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int i = threadIdx.z + blockIdx.z * blockDim.z;

    int Mdim = dim_dev[0];
    int Ndim = dim_dev[1];
    int Kdim = dim_dev[2];

    if ( (i + j + k) % 2 != phase ) return;

    if (i >= inn_dev[0] && i < inn_dev[1] && j >= inn_dev[2] && j < inn_dev[3] && k >= inn_dev[4] && k < inn_dev[5]) {
        u_cur[INDEX(i, j, k)] = 0.5 * ( u_cur[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] - DH2 * f[INDEX(i, j, k)] );
    }

}

template <int phase>
__global__ void rbgs_2D_kernel(type_t * __restrict__ u_cur, type_t * __restrict__ f, const double DH2)
{
    // int k = (threadIdx.x + blockIdx.x * blockDim.x) * 2 + phase + 1;
    int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
    int j = threadIdx.y + blockIdx.y * blockDim.y + 1;
    int i = threadIdx.z + blockIdx.z * blockDim.z;

    int Mdim = dim_dev[0];
    int Ndim = dim_dev[1];
    int Kdim = dim_dev[2];

    
    // if ( j % 2 == 0 ) k = k - 2 * phase + 1;
    if ( (i + j + k) % 2 != phase ) return;

    if (i >= inn_dev[0] && i < inn_dev[1] && j >= inn_dev[2] && j < inn_dev[3] && k >= inn_dev[4] && k < inn_dev[5]) {
        u_cur[INDEX(i, j, k)] = 0.25 * ( u_cur[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] + \
                                         u_cur[INDEX(i, j-1, k)] + u_cur[INDEX(i, j+1, k)] - DH2 * f[INDEX(i, j, k)] );
    }

}

template <int phase>
__global__ void rbgs_3D_kernel(type_t * __restrict__ u_cur, type_t * __restrict__ f, const double DH2)
{
    // int k = (threadIdx.x + blockIdx.x * blockDim.x) * 2 + phase + 1;
    int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
    int j = threadIdx.y + blockIdx.y * blockDim.y + 1;
    int i = threadIdx.z + blockIdx.z * blockDim.z + 1;

    int Mdim = dim_dev[0];
    int Ndim = dim_dev[1];
    int Kdim = dim_dev[2];

    // half of threads; reconsider mesh coloring
    // int offset = 0;
    // if ( j % 2 == 0 ) k = k - 2 * phase + 1;
    // if ( i % 2 == 0 ) {
    //     offset = - 2 * phase + 1;
    //     if ( (k + offset) < inn_dev[4] || (k + offset) >= inn_dev[5])
    //         return;
    // }
    if ( (i + j + k) % 2 != phase ) return;

    if (i >= inn_dev[0] && i < inn_dev[1] && j >= inn_dev[2] && j < inn_dev[3] && k >= inn_dev[4] && k < inn_dev[5]) {
        u_cur[INDEX(i, j, k)] = 1.0 / 6.0 * ( u_cur[INDEX(i, j, k-1)] + u_cur[INDEX(i, j, k+1)] + \
                                              u_cur[INDEX(i, j-1, k)] + u_cur[INDEX(i, j+1, k)] + \
                                              u_cur[INDEX(i-1, j, k)] + u_cur[INDEX(i+1, j, k)] - \
                                              DH2 * f[INDEX(i, j, k)] );
        // u_cur[INDEX(i, j, k) + offset] = 1.0 / 6.0 * ( u_cur[INDEX(i, j, k-1) + offset] + u_cur[INDEX(i, j, k+1) + offset] + \
        //                                                u_cur[INDEX(i, j-1, k) + offset] + u_cur[INDEX(i, j+1, k) + offset] + \
        //                                                u_cur[INDEX(i-1, j, k) + offset] + u_cur[INDEX(i+1, j, k) + offset] - \
        //                                                DH2 * f[INDEX(i, j, k) + offset] );
    }

}

template <int koff, int joff, int ioff>
__global__ void rbgsError_cal_kernel(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ error)
{
    int k = threadIdx.x + blockIdx.x * blockDim.x + koff;
    int j = threadIdx.y + blockIdx.y * blockDim.y + joff;
    int i = threadIdx.z + blockIdx.z * blockDim.z + ioff;

    int Mdim = dim_dev[0];
    int Ndim = dim_dev[1];
    int Kdim = dim_dev[2];

    if (i >= inn_dev[0] && i < inn_dev[1] && j >= inn_dev[2] && j < inn_dev[3] && k >= inn_dev[4] && k < inn_dev[5]) {
        error[INDEX(i, j, k)] = u_new[INDEX(i, j, k)] - u_cur[INDEX(i, j, k)];
    }

}

double iterMethod(type_t *u_cur, type_t *u_new, type_t *error, type_t *f, PARAMS params,
                  INNERS inners, const int rank)
{
    double delta = 0.0;
    int Mdim = params.M;
    int Ndim = params.N;
    int Kdim = params.K;
    int DIM = params.DIM;
    double DH = params.DH;
    double DH2 = DH * DH;
    long long num = Mdim * Ndim * Kdim;

    int istart = inners.istart;     int iend = inners.iend;
    int jstart = inners.jstart;     int jend = inners.jend;
    int kstart = inners.kstart;     int kend = inners.kend;

    // dim3 block(512, 1, 1);
    dim3 block(16, 8, 4);
    dim3 grid;
    grid.x = (kend - kstart + block.x - 1) / block.x;
    grid.y = (jend - jstart + block.y - 1) / block.y;
    grid.z = (iend - istart + block.z - 1) / block.z;

    switch (rank) {
        case 0:
            jacobi_1D_kernel <<<grid, block>>> (u_cur, u_new, error, f, DH2);
            break;
        case 1:
            jacobi_2D_kernel <<<grid, block>>> (u_cur, u_new, error, f, DH2);
            break;
        case 2:
            jacobi_3D_kernel <<<grid, block>>> (u_cur, u_new, error, f, DH2);
            break;
        case 9:
            checkCudaErrors( cudaMemcpy(u_new, u_cur, num * sizeof(type_t), cudaMemcpyDeviceToDevice) );
            // grid.x = ( (kend - kstart) / 2 + 1 + block.x - 1 ) / block.x;    // half of threads
            rbgs_1D_kernel<0> <<<grid, block>>> (u_new, f, DH2);
            rbgs_1D_kernel<1> <<<grid, block>>> (u_new, f, DH2);
            // grid.x = (kend - kstart + block.x - 1) / block.x;
            rbgsError_cal_kernel<1, 0, 0> <<<grid, block>>> (u_cur, u_new, error);
            break;
        case 10:
            checkCudaErrors( cudaMemcpy(u_new, u_cur, num * sizeof(type_t), cudaMemcpyDeviceToDevice) );
            // grid.x = ( (kend - kstart) / 2 + 1 + block.x - 1 ) / block.x;    // half of threads
            rbgs_2D_kernel<0> <<<grid, block>>> (u_new, f, DH2);
            rbgs_2D_kernel<1> <<<grid, block>>> (u_new, f, DH2);
            // grid.x = (kend - kstart + block.x - 1) / block.x;
            rbgsError_cal_kernel<1, 1, 0> <<<grid, block>>> (u_cur, u_new, error);
            break;
        case 11:
            checkCudaErrors( cudaMemcpy(u_new, u_cur, num * sizeof(type_t), cudaMemcpyDeviceToDevice) );
            // grid.x = ( (kend - kstart) / 2 + 1 + block.x - 1 ) / block.x;    // half of threads
            rbgs_3D_kernel<0> <<<grid, block>>> (u_new, f, DH2);
            rbgs_3D_kernel<1> <<<grid, block>>> (u_new, f, DH2);
            // grid.x = (kend - kstart + block.x - 1) / block.x;
            rbgsError_cal_kernel<1, 1, 1> <<<grid, block>>> (u_cur, u_new, error);
            break;

    }


    checkCudaErrors( cublasDnrm2(handle, num, error, 1, &delta) );
    return delta;

}
