#pragma once

#include "helper_cuda.h"
#include "helper_string.h"
#include <iostream>

class DisableCopy
{
public:
    DisableCopy() = default;
    DisableCopy(DisableCopy const &) = delete;
    DisableCopy &operator=(DisableCopy const &) = delete;
};


template <typename T>
class CudaArray: public DisableCopy
{
private:
    cudaArray_t m_cuArray{};
    uint3 m_dim{};

public:
    explicit CudaArray(uint3 const &_dim)
        : m_dim(_dim)
    {
        std::cout << m_dim.x << std::endl;
        cudaExtent extent = make_cudaExtent(m_dim.x, m_dim.y, m_dim.z);
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<T>();
        checkCudaErrors(cudaMalloc3DArray(&m_cuArray, &channelDesc, extent, cudaArraySurfaceLoadStore));
    }

    ~CudaArray()
    {
        checkCudaErrors(cudaFreeArray(m_cuArray));
    }

    void copyToDevice(T const *_data);

    void copyToHost(T const *_data);

    cudaArray_t getArray() const
    {
        return m_cuArray;
    }

};

template <typename T>
void CudaArray<T>::copyToDevice(T const *_data)
{
    std::cout << "copyToDevice" << std::endl;
    cudaMemcpy3DParms copy3DParams = {0};
    copy3DParams.srcPtr = make_cudaPitchedPtr((void *)_data, sizeof(T) * m_dim.x, m_dim.x, m_dim.y);
    copy3DParams.dstArray = m_cuArray;
    copy3DParams.extent = make_cudaExtent(m_dim.x, m_dim.y, m_dim.z);
    copy3DParams.kind = cudaMemcpyHostToDevice;
    checkCudaErrors(cudaMemcpy3D(&copy3DParams));
}

template <typename T>
void CudaArray<T>::copyToHost(T const *_data)
{
    std::cout << "copyToHost" << std::endl;
    cudaMemcpy3DParms copy3DParams = {};
    copy3DParams.srcArray = m_cuArray;
    copy3DParams.dstPtr = make_cudaPitchedPtr((void *)_data, sizeof(T) * m_dim.x, m_dim.x, m_dim.y);
    copy3DParams.extent = make_cudaExtent(m_dim.x, m_dim.y, m_dim.z);
    copy3DParams.kind = cudaMemcpyDeviceToHost;
    checkCudaErrors(cudaMemcpy3D(&copy3DParams));
}


// Accessor pattern
template <typename T>
class CudaSurfaceAccessor
{
private:
    cudaSurfaceObject_t m_cuSurf{};

public:
    explicit CudaSurfaceAccessor(cudaSurfaceObject_t m_cuSurf)
        : m_cuSurf(m_cuSurf)
    {}

    template <cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap>
    __device__ __forceinline__ T read(int x, int y, int z) const
    {
        return surf3Dread<T>(m_cuSurf, x * sizeof(T), y, z, mode);  // surface object use byte coordinare x
    }

    template <cudaSurfaceBoundaryMode mode = cudaBoundaryModeTrap>
    __device__ __forceinline__ void write(T val, int x, int y, int z) const
    {
        surf3Dwrite<T>(val, m_cuSurf, x * sizeof(T), y, z, mode);
    }

};


template <typename T>
class CudaSurface: public CudaArray<T>
{
private:
    cudaSurfaceObject_t m_cuSurf{};

public:
    explicit CudaSurface(uint3 const &_dim)
        : CudaArray<T>(_dim)
    {
        cudaResourceDesc resDesc{};
        resDesc.resType = cudaResourceTypeArray;
        resDesc.res.array.array = CudaArray<T>::getArray();
        checkCudaErrors(cudaCreateSurfaceObject(&m_cuSurf, &resDesc));
    }

    ~CudaSurface()
    {
        checkCudaErrors(cudaDestroySurfaceObject(m_cuSurf));
    }

    cudaSurfaceObject_t getSurface() const
    {
        return m_cuSurf;
    }
    
    CudaSurfaceAccessor<T> accessSurface() const
    {
        return CudaSurfaceAccessor<T>{m_cuSurf};
    }

};








template <typename T>
class CudaTextureAccessor
{
private:
    cudaTextureObject_t m_cuTex;

public:
    __device__ __forceinline__ T sample(float x, float y, float z)
    {
        return tex3D<T>(m_cuTex, x, y, z);
    }

};


// surface object and texture object can be bound to the same underlying cudaArray
template <typename T>
class CudaTexture: public CudaSurface<T>
{
private:
    struct Parameters {
        cudaTextureAddressMode addressMode{cudaAddressModeBorder};
        cudaTextureFilterMode filterMode{cudaFilterModeLinear};
        cudaTextureReadMode readMode{cudaReadModeElementType};
        int normalizedCoords{false};
    };

    cudaTextureObject_t m_cuTex{};

public:
    explicit CudaTexture(uint3 const &_dim, Parameters const &_args = {})
        : CudaSurface<T>(_dim)
    {
        cudaResourceDesc resDesc{};
        resDesc.resType = cudaResourceTypeArray;
        resDesc.res.array.array = CudaSurface<T>::getArray();

        cudaTextureDesc texDesc{};
        texDesc.addressMode[0] = _args.addressMode;
        texDesc.addressMode[1] = _args.addressMode;
        texDesc.addressMode[2] = _args.addressMode;
        texDesc.filterMode = _args.filterMode;
        texDesc.readMode = _args.readMode;
        texDesc.normalizedCoords = _args.normalizedCoords;

        checkCudaErrors(cudaCreateTextureObject(&m_cuTex, &resDesc, &texDesc, NULL));
    }

    ~CudaTexture()
    {
        checkCudaErrors(cudaDestroyTextureObject(m_cuTex));
    }

    cudaTextureObject_t getTexture() const
    {
        return m_cuTex;
    }

    CudaTextureAccessor<T> accessTexture() const
    {
        return {m_cuTex};
    }

};




