// #pragma once

#include "CudaArray.cuh"

#include <cmath>
#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <chrono>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "helper_cuda.h"
#include "helper_string.h"
#include <cuda_runtime.h>
#include <cublas_v2.h>

const double PI = 3.14159265535;

/*
 * ix  -> the index
 * nx  -> the dimension
 * ix0 -> the offset from beginning
 * ixn -> the offset from ending
 */
__host__ __device__ __forceinline__
bool inbounds(const int ix, const int iy, const int iz,
              const int nx, const int ny, const int nz,
              const int ix0 = 0, const int iy0 = 0, const int iz0 = 0,
              const int ixn = 0, const int iyn = 0, const int izn = 0)
{
    if (ix < ix0 || iy < iy0 || iz < iz0) return false;
    if (ix > nx - 1 - ixn || iy > ny - 1 - iyn || iz > nz - 1 - izn) return false;
    return true;
}

template <typename T>
struct Functor
{
    __host__ __device__ T operator()(int x, int y, int z, T dh, int fhz) const
    {
        return (- 3.0 * fhz * fhz * PI * PI * std::sin(fhz * PI * x * dh)
                                            * std::sin(fhz * PI * y * dh)
                                            * std::sin(fhz * PI * z * dh));
    }

};


template <typename T>
__global__ void initial_kernel(CudaSurfaceAccessor<T> u, CudaSurfaceAccessor<T> f,
                               const int Xdim, const int Ydim, const int Zdim,
                               T const dh, Functor<T> functor, const int fhz)
{
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    int z = threadIdx.z + blockIdx.z * blockDim.z;

    if (!inbounds(x, y, z, Xdim, Ydim, Zdim)) return;

    T uVal = 0.0;
    T fVal = functor(x, y, z, dh, fhz);
    // printf("%f\t", fVal);
    u.write(uVal, x, y, z);
    f.write(fVal, x, y, z);
}

// template <typename T>
// __global__ void calcuDelta_kernel(CudaSurfaceAccessor<T> u, CudaSurfaceAccessor<T> f, CudaSurfaceAccessor<T> res,
//                                   const int nx, const int ny, const int nz, T const ddhh)
// {
//     int x = threadIdx.x + blockIdx.x * blockDim.x + 1;
//     int y = threadIdx.y + blockIdx.y * blockDim.y + 1;
//     int z = threadIdx.z + blockIdx.z * blockDim.z + 1;

//     if (!inbounds(x, y, z, nx, ny, nz, 1, 1, 1, 1, 1, 1)) return;

//     // p-previous   c-current   n-next   t-next time
//     T uc = u.read<cudaBoundaryModeClamp>(x, y, z) * 6;
//     T xp = u.read<cudaBoundaryModeClamp>(x - 1, y, z);
//     T xn = u.read<cudaBoundaryModeClamp>(x + 1, y, z);
//     T yp = u.read<cudaBoundaryModeClamp>(x, y - 1, z);
//     T yn = u.read<cudaBoundaryModeClamp>(x, y + 1, z);
//     T zp = u.read<cudaBoundaryModeClamp>(x, y, z - 1);
//     T zn = u.read<cudaBoundaryModeClamp>(x, y, z + 1);

//     T fc = f.read(x, y, z);
//     T resVal = fc - (xp + xn + yp + yn + zp + zn - uc) / ddhh;
//     res.write(resVal, x, y, z);

// }


template <typename T, int phase>
__global__ void rbgs_kernel(CudaSurfaceAccessor<T> u, CudaSurfaceAccessor<T> f,
                            const int nx, const int ny, const int nz, T const ddhh)
{
    int x = threadIdx.x + blockIdx.x * blockDim.x + 1;
    int y = threadIdx.y + blockIdx.y * blockDim.y + 1;
    int z = threadIdx.z + blockIdx.z * blockDim.z + 1;

    if (!inbounds(x, y, z, nx, ny, nz, 1, 1, 1, 1, 1, 1)) return;
    if ((x + y + z) % 2 != phase) return;

    // p-previous   c-current   n-next   t-next time
    T xp = u.read(x - 1, y, z);
    T xn = u.read(x + 1, y, z);
    T yp = u.read(x, y - 1, z);
    T yn = u.read(x, y + 1, z);
    T zp = u.read(x, y, z - 1);
    T zn = u.read(x, y, z + 1);

    T fc = f.read(x, y, z) * ddhh;
    T ut = (xp + xn + yp + yn + zp + zn - fc) * static_cast<T>(1.0 / 6.0);
    // printf("%f\t", ut); 
    u.write(ut, x, y, z);

}





template <typename T>
class PossionSolver: public DisableCopy
{
private:
    // inline static dim3 block{32, 16, 1};
    
    unsigned int Xdim;  // fastest axis
    unsigned int Ydim;  // subfast axis
    unsigned int Zdim;  // slowest axis
    unsigned int n_levels;
    unsigned int n_cycles;
    T DH;
    unsigned int fhz;
    Functor<T> functor;
    inline static dim3 block{32, 16, 1};    // c++17
    std::unique_ptr<CudaSurface<T>> u;
    std::unique_ptr<CudaSurface<T>> f;
    std::unique_ptr<CudaSurface<T>> delta_res;
    std::vector<std::unique_ptr<CudaSurface<T>>> res;
    std::vector<std::unique_ptr<CudaSurface<T>>> res2;
    std::vector<std::unique_ptr<CudaSurface<T>>> err2;
    std::vector<uint3> dimInfos;
    std::vector<T> dhInfos;
    std::vector<T> deltaVector{};

public:
    explicit PossionSolver(unsigned int K = 128, unsigned int N = 128, unsigned int M = 128,
                           unsigned int levels = 4, unsigned int cycles = 16,
                           T dh = static_cast<T>(1.0 / 127), int hz = 4)
        : Xdim(K), Ydim(N), Zdim(M), n_levels(levels), n_cycles(cycles), DH(dh), fhz(hz)
        , functor{}
        , u(std::make_unique<CudaSurface<T>>(uint3{Xdim, Ydim, Zdim}))  // c++14
        , f(std::make_unique<CudaSurface<T>>(uint3{Xdim, Ydim, Zdim}))
        , delta_res(std::make_unique<CudaSurface<T>>(uint3{Xdim, Ydim, Zdim}))
    {
        unsigned int x_lev = Xdim;
        unsigned int y_lev = Ydim;
        unsigned int z_lev = Zdim;
        T dh_lev = DH;
        dimInfos.emplace_back(uint3{x_lev, y_lev, z_lev});
        dhInfos.emplace_back(dh_lev);
        for (int i = 0; i < n_levels - 1; i++) {
            res.emplace_back(std::make_unique<CudaSurface<T>>(uint3{x_lev, y_lev, z_lev}));
            x_lev /= 2;
            y_lev /= 2;
            z_lev /= 2;
            dh_lev *= 2;

            res2.emplace_back(std::make_unique<CudaSurface<T>>(uint3{x_lev, y_lev, z_lev}));
            err2.emplace_back(std::make_unique<CudaSurface<T>>(uint3{x_lev, y_lev, z_lev}));
            dimInfos.emplace_back(uint3{x_lev, y_lev, z_lev});
            dhInfos.emplace_back(dh_lev);
        }

    }

    // initialization on device
    void initialization()
    {
        // dim3 block{32, 16, 1};
        dim3 grid;
        grid.x = (Xdim + block.x - 1) / block.x;
        grid.y = (Ydim + block.y - 1) / block.y;
        grid.z = (Zdim + block.z - 1) / block.z;

        initial_kernel<T><<<grid, block>>>(u->accessSurface(), f->accessSurface(), Xdim, Ydim, Zdim, DH, functor, fhz);
        checkCudaErrors(cudaGetLastError());
    }

    // initialization on host and memcpy to device
    void initialization(T const *u_host, T const *f_host)
    {
        u->copyToDevice(u_host);
        f->copyToDevice(f_host);
    }

    // output device data to host memory
    void outputData(T const *u_host, T const *f_host)
    {
        u->copyToHost(u_host);
        f->copyToHost(f_host);

        FILE *fp;
        char filename[128];
        sprintf(filename, "../output/MultigridCPP_%d_%dD_GPU.bin", n_levels, 3);
        fp = fopen(filename, "wb");
        fwrite(u_host, sizeof(T), Xdim*Ydim*Zdim, fp);
        fclose(fp);

        FILE *fp1;
        char filename1[128];
        sprintf(filename1, "../output/MultigridCPP_%d_%dD_dim_type_GPU.txt", n_levels, 3);
        fp1 = fopen(filename1, "w");
        fprintf(fp1, "DIM %d %d %d sizeof %d", Xdim, Ydim, Zdim, sizeof(T));
        fclose(fp1);

        FILE *fp2;
        char filename2[128];
        sprintf(filename2, "../output/MultigridCPP_%d_%dD_delta_GPU.bin", n_levels, 3);
        fp2 = fopen(filename2, "wb");
        fwrite(deltaVector.data(), sizeof(T), deltaVector.size(), fp2);
        fclose(fp2);

    }

    // smoother
    void smooth(CudaSurfaceAccessor<T> u, CudaSurfaceAccessor<T> f, const int level);

    void vcycle(CudaSurfaceAccessor<T> u, CudaSurfaceAccessor<T> f, const int level)
    {
        if (level >= n_levels - 1) {
            smooth(u, f, level);
            return;
        }

        // checkCudaErrors(cudaDeviceSynchronize());


    }

    void multigrid()
    {
        cublasHandle_t handle;
        cublasCreate(&handle);
        
        int iter = 0;
        T eps = 1.0e-6;
        T delta = 1.0;

        using double_ms = std::chrono::duration<double, std::milli>;
        auto time_start = std::chrono::steady_clock::now();
        dim3 grid;
        grid.x = (Xdim - 2 + block.x - 1) / block.x;
        grid.y = (Ydim - 2 + block.y - 1) / block.y;
        grid.z = (Zdim - 2 + block.z - 1) / block.z;
        // calcuDelta_kernel<<<grid, block>>>
        //     (u.get()->accessSurface(), f.get()->accessSurface(), delta_res.get()->accessSurface(), Xdim, Ydim, Zdim, DH * DH);

        // if constexpr (std::is_same_v<T, float>) // c++17
        //     checkCudaErrors(cublasSnrm2(handle, Xdim*Ydim*Zdim, delta_res.get()->getArray(), 1, &delta));
        // else if constexpr (std::is_same_v<T, double>)
        //     checkCudaErrors(cublasDnrm2(handle, Xdim*Ydim*Zdim, delta_res.get()->getArray(), 1, &delta));
        // else
        //     std::cout << "Note: data type are neither float or double!!!" << std::endl;
        
        deltaVector.emplace_back(delta);
        auto time_inter = std::chrono::steady_clock::now();
        printf("iter: %6d\tdelta=%15g\ttime:%fms\n", iter, delta,
                std::chrono::duration_cast<double_ms>(time_inter - time_start).count());

        for (int loop = 0; loop < 20; loop++) {
        // while (delta > eps) {
            vcycle(u.get()->accessSurface(), f.get()->accessSurface(), 0);

            // calcuDelta_kernel<<<grid, block>>>
            //     (u.get()->accessSurface(), f.get()->accessSurface(), delta_res.get()->accessSurface(), Xdim, Ydim, Zdim, DH * DH);
            // if constexpr (std::is_same_v<T, float>) // c++17
            //     checkCudaErrors(cublasSnrm2(handle, Xdim*Ydim*Zdim, delta_res.get()->getArray(), 1, &delta));
            // else if constexpr (std::is_same_v<T, double>)
            //     checkCudaErrors(cublasDnrm2(handle, Xdim*Ydim*Zdim, delta_res.get()->getArray(), 1, &delta));
            // else
            //     std::cout << "Note: data type are neither float or double!!!" << std::endl;

            deltaVector.emplace_back(delta);
            iter += (n_cycles) * (n_levels > 1 ? (2 * n_levels - 1) : 1);

            time_inter = std::chrono::steady_clock::now();
            printf("iter: %6d\tdelta=%15g\ttime:%fms\n", iter, delta,
                    std::chrono::duration_cast<double_ms>(time_inter - time_start).count());


        }

        auto time_end = std::chrono::steady_clock::now();
        printf("\nElapsed time: %fms\n",
                std::chrono::duration_cast<double_ms>(time_end - time_start).count());

        return;
    }


};


template <typename T>
void PossionSolver<T>::smooth(CudaSurfaceAccessor<T> u, CudaSurfaceAccessor<T> f, const int level)
{
    int nx = dimInfos[level].x;
    int ny = dimInfos[level].y;
    int nz = dimInfos[level].z;
    T ddhh = dhInfos[level] * dhInfos[level];

    dim3 grid;
    grid.x = (nx - 2 + block.x - 1) / block.x;
    grid.y = (ny - 2 + block.y - 1) / block.y;
    grid.z = (nz - 2 + block.z - 1) / block.z;

    for (int iter = 0; iter < n_cycles; iter++) {
        rbgs_kernel<T, 0><<<grid, block>>>(u, f, nx, ny, nz, ddhh);
        rbgs_kernel<T, 1><<<grid, block>>>(u, f, nx, ny, nz, ddhh);
    }

}




int main(int argc, char **argv)
{
    const unsigned int n = 400;
    PossionSolver<float> myTest(n, n, n, 1, 16, 1.0/(n - 1), 4);

    float *u_host = new float [n * n * n];
    float *f_host = new float [n * n * n];

    // myTest.initialization(u_host, f_host);
    myTest.initialization();

    myTest.multigrid();

    myTest.outputData(u_host, f_host);

    // for (int i = 0; i < n * n * n; i++) {
    //     std::cout << f_host[i] << std::endl;
    // }


    return 0;
}


