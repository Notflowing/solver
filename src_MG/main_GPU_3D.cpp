#include "header.h"

const char *iterMth[] = {"DEFAULT", "JACOBI", "GAUSSSEIDEL", "SOR", "RBGS", "MGRBGS"};

// const double relax_rbgs = 0.5;
const double relax_rbgs = 1.0;


// Set the number of levels and the number of cycles
const int n_levels = 4; // 0, 1, 2, 3--total 4; level 3 is the base solver
const int n_cycles = 16;
type_t *res [n_levels - 1];
type_t *res2[n_levels - 1];
type_t *err2[n_levels - 1];
DIMINFO dimsinfo[n_levels];
INNERS innersinfo[n_levels];

__constant__ DIMINFO dimsinfo_dev[n_levels];
__constant__ INNERS innersinfo_dev[n_levels];


typedef void (*MGIterSolver_t)(type_t *u, type_t *f, const DIMINFO diminfo, const INNERS inners);
MGIterSolver_t MGIterSolver;


void paramsCheck(int argc, char **argv, int *M, int *N, int *K, int *DIM,
                 int *method, const char **iterMth)
{
    switch (argc) {
        case 1:
            *M = *N = 1;
            *K = 1024;
            *DIM = 1;
            *method = 1;
            printf("Default dimension: 1024 and default method: JACOBI\n");
            break;
        case 2:
            *method = atoi(argv[1]);
            *M = *N = 1;
            *K = 1024;
            *DIM = 1;
            printf("1 dimension: %d and iterative method:%s\n", *K, iterMth[*method]);
            break;
        case 3:
            *method = atoi(argv[1]);
            *M = *N = 1;
            *K = atoi(argv[2]);
            *DIM = 1;
            printf("1 dimension: %d and iterative method:%s\n", *K, iterMth[*method]);
            break;
        case 4:
            *method = atoi(argv[1]);
            *M = 1;
            *N = atoi(argv[2]);
            *K = atoi(argv[3]);
            *DIM = 2;
            printf("2 dimension: %d x %d and iterative method:%s\n", *N, *K, iterMth[*method]);
            break;
        case 5:
            *method = atoi(argv[1]);
            *M = atoi(argv[2]);
            *N = atoi(argv[3]);
            *K = atoi(argv[4]);
            *DIM = 3;
            printf("3 dimension: %d x %d x %d and iterative method:%s\n", *M, *N, *K, iterMth[*method]);
            break;
        default:
            printf("Error! Usage: ./bin/solver [method] [M] [N] [K]!!!\n");
            exit(1);
    }

}


void initialize(type_t *u, type_t *f, const int Mdim, const int Ndim, const int Kdim,
                const int DIM, const double DH)
{
    int i, j, k;
    long long idx;
    // auto func = [](const int i, const int j, const int k, const double DH){return sin(KHZ * PI * k * DH);};
    std::function<type_t(const int, const int, const int, const double)> func;
    switch (DIM) {
        case (1):
            func = [](const int i, const int j, const int k, const double DH) \
                   {return  - 1 * KHZ * KHZ * PI * PI * sin(KHZ * PI * k * DH);};
            break;
        case (2):
            func = [](const int i, const int j, const int k, const double DH) \
                   {return  - 2 * KHZ * KHZ * PI * PI * sin(KHZ * PI * k * DH) * sin(KHZ * PI * j * DH);};
            break;
        case (3):
            func = [](const int i, const int j, const int k, const double DH) \
                   {return  - 3 * KHZ * KHZ * PI * PI * sin(KHZ * PI * k * DH) * sin(KHZ * PI * j * DH) * sin(KHZ * PI * i * DH);};
            break;
    }
    
    for (i = 0; i < Mdim; i++) {
        for (j = 0; j < Ndim; j++) {
            for (k = 0; k < Kdim; k++) {
                idx = INDEX(i, j, k);
                u[idx] = 0.0;
                f[idx] = func(i, j, k, DH);
                // printf("%f\t", f[idx]);
            }
        }
    }
    
}

void writeResult(type_t *u, std::vector<double> &delta,
                 const int M, const int N, const int K, const int DIM, double totalTime)
{

    type_t *u_host = (type_t *)malloc(sizeof(type_t) * M*N*K);
    checkCudaErrors( cudaMemcpy(u_host, u, sizeof(type_t) * M*N*K, cudaMemcpyDeviceToHost) );

    FILE *fp;
    char filename[128];
    sprintf(filename, "output/Multigrid_%d_%dD_GPU.bin", n_levels, DIM);
    fp = fopen(filename, "wb");
    fwrite(u_host, sizeof(type_t), M*N*K, fp);
    fclose(fp);

    FILE *fp1;
    char filename1[128];
    sprintf(filename1, "output/Multigrid_%d_%dD_dim_type_GPU.txt", n_levels, DIM);
    fp1 = fopen(filename1, "w");
    fprintf(fp1, "DIM %d %d %d sizeof %d totalTime %lf", M, N, K, sizeof(type_t), totalTime);
    fclose(fp1);

    FILE *fp2;
    char filename2[128];
    sprintf(filename2, "output/Multigrid_%d_%dD_delta_GPU.bin", n_levels, DIM);
    fp2 = fopen(filename2, "wb");
    fwrite(delta.data(), sizeof(double), delta.size(), fp2);
    fclose(fp2);

}


inline type_t *allocMem_dev(const int M, const int N, const int K) {
    size_t num = M * N * K;
    type_t *ptr;
    checkCudaErrors( cudaMalloc((void **)&ptr, sizeof(type_t) * num) );
    checkCudaErrors( cudaMemset(ptr, 0, sizeof(type_t) * num) );

    return ptr;
}



template <int phase>
__global__ void rbgs_kernel(type_t * __restrict__ u, type_t * __restrict__ f, const int level)
{
    int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
    int j = threadIdx.y + blockIdx.y * blockDim.y + 1;
    int i = threadIdx.z + blockIdx.z * blockDim.z + 1;

    DIMINFO diminfo = dimsinfo_dev[level];
    INNERS inners = innersinfo_dev[level];

    int Mdim = diminfo.M;
    int Ndim = diminfo.N;
    int Kdim = diminfo.K;
    double DH2 = diminfo.DH * diminfo.DH;

    // printf("%d_%d_%d\n", Mdim, Ndim, Kdim);
    // printf("%d_%d_%d_%d_%d_%d\n%d_%d_%d_%e\n", inners.istart, inners.iend,
    //                                            inners.jstart, inners.jend,
    //                                            inners.kstart, inners.kend,
    //                                            diminfo.M, diminfo.N, diminfo.K, diminfo.DH);

    if ( (i + j + k) % 2 != phase ) return;

    // =====================================
    // !!!!!! must add a relax factor !!!!!!
    // =====================================
    if (i >= inners.istart && i < inners.iend && j >= inners.jstart && j < inners.jend && k >= inners.kstart && k < inners.kend) {
        u[INDEX(i, j, k)] = relax_rbgs * (1.0 / 6.0) * ( u[INDEX(i, j, k-1)] + u[INDEX(i, j, k+1)] + \
                                                         u[INDEX(i, j-1, k)] + u[INDEX(i, j+1, k)] + \
                                                         u[INDEX(i-1, j, k)] + u[INDEX(i+1, j, k)] - \
                                                         DH2 * f[INDEX(i, j, k)] ) + (1 - relax_rbgs) * u[INDEX(i, j, k)];
    }

}

void smooth(type_t *u, type_t *f, const int level)
{
    dim3 block(32, 16, 1);
    dim3 grid;
    grid.x = (innersinfo[level].kend - innersinfo[level].kstart + block.x - 1) / block.x;
    grid.y = (innersinfo[level].jend - innersinfo[level].jstart + block.y - 1) / block.y;
    grid.z = (innersinfo[level].iend - innersinfo[level].istart + block.z - 1) / block.z;
    int iter;
    for (iter = 0; iter < n_cycles; iter++) {
        rbgs_kernel<0><<<grid, block>>>(u, f, level);
        rbgs_kernel<1><<<grid, block>>>(u, f, level);
    }

}


void MGgaussseidel_2D(type_t *u, type_t *f, const DIMINFO diminfo, const INNERS inners)
{
    int Mdim = diminfo.M;
    int Ndim = diminfo.N;
    int Kdim = diminfo.K;
    double DH2 = diminfo.DH * diminfo.DH;
// printf("%s, %d\n", __FILE__, __LINE__);
    int i, j, k;
    for (i = inners.istart; i < inners.iend; i++) {
        for (j = inners.jstart; j < inners.jend; j++) {
            for (k = inners.kstart; k < inners.kend; k++) {
                u[INDEX(i, j, k)] = 0.25 * ( u[INDEX(i, j, k-1)] + u[INDEX(i, j, k+1)] + \
                                             u[INDEX(i, j-1, k)] + u[INDEX(i, j+1, k)] - \
                                             DH2 * f[INDEX(i, j, k)] );
            }
        }
    }

}



__global__ void residual_kernel(type_t * __restrict__ r, type_t * __restrict__ u, type_t * __restrict__ f, const int level)
{
    int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
    int j = threadIdx.y + blockIdx.y * blockDim.y + 1;
    int i = threadIdx.z + blockIdx.z * blockDim.z + 1;

    DIMINFO diminfo = dimsinfo_dev[level];
    INNERS inners = innersinfo_dev[level];

    int Mdim = diminfo.M;
    int Ndim = diminfo.N;
    int Kdim = diminfo.K;
    double DH2 = diminfo.DH * diminfo.DH;

    if (i >= inners.istart && i < inners.iend && j >= inners.jstart && j < inners.jend && k >= inners.kstart && k < inners.kend) {
        r[INDEX(i, j, k)] = f[INDEX(i, j, k)] - ( u[INDEX(i, j, k-1)] + u[INDEX(i, j, k+1)] + \
                                                  u[INDEX(i, j-1, k)] + u[INDEX(i, j+1, k)] + \
                                                  u[INDEX(i-1, j, k)] + u[INDEX(i+1, j, k)] - \
                                                  6 * u[INDEX(i, j, k)] ) / DH2;
    }

}



__global__ void restriction_kernel(type_t * __restrict__ r2, type_t * __restrict__ r, type_t * __restrict__ e2, const int level)
{
    int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
    int j = threadIdx.y + blockIdx.y * blockDim.y + 1;
    int i = threadIdx.z + blockIdx.z * blockDim.z + 1;

    DIMINFO diminfo_coar = dimsinfo_dev[level];
    DIMINFO diminfo_fine = dimsinfo_dev[level - 1];

    int Mdim_coar = diminfo_coar.M;
    int Ndim_coar = diminfo_coar.N;
    int Kdim_coar = diminfo_coar.K;

    int Mdim_fine = diminfo_fine.M;
    int Ndim_fine = diminfo_fine.N;
    int Kdim_fine = diminfo_fine.K;

    // if (Ndim_coar != Kdim_coar || Ndim_fine != Kdim_fine || Ndim_coar * 2 != Ndim_fine) {
    //     printf("restrict\n%d_%d_%d_%d_%d_%d\n", Mdim_coar, Ndim_coar, Kdim_coar, Mdim_fine, Ndim_fine, Kdim_fine);
    // }

    if (i >= 1 && i < Mdim_coar - 1 && j >= 1 && j < Ndim_coar - 1 && k >= 1 && k < Kdim_coar - 1) {
        r2[INDEX_coar(i, j, k)] = (1.0 / 12.0) * (6 * r[INDEX_fine(2*i, 2*j, 2*k)] + \
                                                      r[INDEX_fine(2*i, 2*j, 2*k-1)] + r[INDEX_fine(2*i, 2*j, 2*k+1)] + \
                                                      r[INDEX_fine(2*i, 2*j-1, 2*k)] + r[INDEX_fine(2*i, 2*j+1, 2*k)] + \
                                                      r[INDEX_fine(2*i-1, 2*j, 2*k)] + r[INDEX_fine(2*i+1, 2*j, 2*k)]);
        e2[INDEX_coar(i, j, k)] = 0.0;
    }

}



__global__ void prolongate_kernel(type_t * __restrict__ u, type_t * __restrict__ e2, const int level)
{
    int k = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int i = threadIdx.z + blockIdx.z * blockDim.z;

    DIMINFO diminfo_coar = dimsinfo_dev[level];
    DIMINFO diminfo_fine = dimsinfo_dev[level - 1];

    int Mdim_coar = diminfo_coar.M;
    int Ndim_coar = diminfo_coar.N;
    int Kdim_coar = diminfo_coar.K;

    int Mdim_fine = diminfo_fine.M;
    int Ndim_fine = diminfo_fine.N;
    int Kdim_fine = diminfo_fine.K;

    // if (Ndim_coar != Kdim_coar || Ndim_fine != Kdim_fine || Ndim_coar * 2 != Ndim_fine) {
    //     printf("prolonge\n%d_%d_%d_%d_%d_%d\n", Mdim_coar, Ndim_coar, Kdim_coar, Mdim_fine, Ndim_fine, Kdim_fine);
    // }
    // printf("%d_%d_%d_%d_%d_%d\n", Mdim_coar, Ndim_coar, Kdim_coar, Mdim_fine, Ndim_fine, Kdim_fine);


    if (i >= 0 && i < (Mdim_coar - 1) && j >= 0 && j < (Ndim_coar - 1) && k >= 0 && k < (Kdim_coar - 1)) {
        u[INDEX_fine(2*i, 2*j, 2*k)]     += e2[INDEX_coar(i, j, k)];

        u[INDEX_fine(2*i, 2*j, 2*k+1)]   += 0.5 * (e2[INDEX_coar(i, j, k)] + e2[INDEX_coar(i, j, k+1)]);
        u[INDEX_fine(2*i, 2*j+1, 2*k)]   += 0.5 * (e2[INDEX_coar(i, j, k)] + e2[INDEX_coar(i, j+1, k)]);
        u[INDEX_fine(2*i+1, 2*j, 2*k)]   += 0.5 * (e2[INDEX_coar(i, j, k)] + e2[INDEX_coar(i+1, j, k)]);

        u[INDEX_fine(2*i, 2*j+1, 2*k+1)] += 0.25* (e2[INDEX_coar(i, j, k)] + e2[INDEX_coar(i, j, k+1)] + \
                                                   e2[INDEX_coar(i, j+1, k)] + e2[INDEX_coar(i, j+1, k+1)]);
        u[INDEX_fine(2*i+1, 2*j, 2*k+1)] += 0.25* (e2[INDEX_coar(i, j, k)] + e2[INDEX_coar(i, j, k+1)] + \
                                                   e2[INDEX_coar(i+1, j, k)] + e2[INDEX_coar(i+1, j, k+1)]);
        u[INDEX_fine(2*i+1, 2*j+1, 2*k)] += 0.25* (e2[INDEX_coar(i, j, k)] + e2[INDEX_coar(i, j+1, k)] + \
                                                   e2[INDEX_coar(i+1, j, k)] + e2[INDEX_coar(i+1, j+1, k)]);

        u[INDEX_fine(2*i+1, 2*j+1, 2*k+1)] += 0.125* (e2[INDEX_coar(i, j, k)] + e2[INDEX_coar(i, j, k+1)] + \
                                                      e2[INDEX_coar(i, j+1, k)] + e2[INDEX_coar(i, j+1, k+1)] + \
                                                      e2[INDEX_coar(i+1, j, k)] + e2[INDEX_coar(i+1, j, k+1)] + \
                                                      e2[INDEX_coar(i+1, j+1, k)] + e2[INDEX_coar(i+1, j+1, k+1)]);
    }

}



__global__ void calcuDelta_dev(type_t *u, type_t *f, type_t *res, const int level)
{
    DIMINFO diminfo = dimsinfo_dev[level];
    INNERS inners = innersinfo_dev[level];

    int Mdim = diminfo.M;
    int Ndim = diminfo.N;
    int Kdim = diminfo.K;
    double DH2 = diminfo.DH * diminfo.DH;

    int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
    int j = threadIdx.y + blockIdx.y * blockDim.y + 1;
    int i = threadIdx.z + blockIdx.z * blockDim.z + 1;

    if (i >= inners.istart && i < inners.iend && j >= inners.jstart && j < inners.jend && k >= inners.kstart && k < inners.kend) {
        res[INDEX(i, j, k)] = f[INDEX(i, j, k)] - ( u[INDEX(i, j, k-1)] + u[INDEX(i, j, k+1)] + \
                                                    u[INDEX(i, j-1, k)] + u[INDEX(i, j+1, k)] + \
                                                    u[INDEX(i-1, j, k)] + u[INDEX(i+1, j, k)] - \
                                                    6 * u[INDEX(i, j, k)] ) / DH2;
    }

}

void vcycle(type_t *u, type_t *f, const int level)
{
    if (level >= n_levels - 1) {

        smooth(u, f, level);

        // smooth_gs(u, f, level);
        return;
    }

    dim3 block(32, 16, 1);
    dim3 grid;
    dim3 grid_1;
    dim3 grid_2;
    grid.x = (innersinfo[level].kend - innersinfo[level].kstart + block.x - 1) / block.x;
    grid.y = (innersinfo[level].jend - innersinfo[level].jstart + block.y - 1) / block.y;
    grid.z = (innersinfo[level].iend - innersinfo[level].istart + block.z - 1) / block.z;

    type_t *r = res[level];
    type_t *r2 = res2[level];
    type_t *e2 = err2[level];
    // Pre-smoothing
    smooth(u, f, level);

    // Compute residual error
    residual_kernel <<<grid, block>>> (r, u, f, level);


    // Restriction and fill zero
    grid_1.x = (dimsinfo[level + 1].K - 2 + block.x - 1) / block.x;
    grid_1.y = (dimsinfo[level + 1].N - 2 + block.y - 1) / block.y;
    grid_1.z = (dimsinfo[level + 1].M - 2 + block.z - 1) / block.z;
    restriction_kernel <<<grid_1, block>>> (r2, r, e2, level + 1);


    // Recursion vcycle
    vcycle(e2, r2, level + 1);


    // Prolongation and correction
    grid_2.x = (dimsinfo[level + 1].K - 1 + block.x - 1) / block.x;
    grid_2.y = (dimsinfo[level + 1].N - 1 + block.y - 1) / block.y;
    grid_2.z = (dimsinfo[level + 1].M - 1 + block.z - 1) / block.z;
    prolongate_kernel <<<grid_2, block>>> (u, e2, level + 1);


    // Post-smoothing
    smooth(u, f, level);
    // checkCudaErrors( cudaDeviceSynchronize() );

    return;

}


void multigrid(type_t *u, type_t *f,
               const int M, const int N, const int K,
               const int DIM, const double DH)
{

    int i;
    dimsinfo[0] = {M, N, K, DH};
    int m = M;
    int n = N;
    int k = K;
    double dh = DH;
    // allocate memory of residual ans error for different level
    for (i = 0; i < n_levels - 1; i++) {

        res[i]  = allocMem_dev(m, n, k);
        m /= 2;
        n /= 2;
        k /= 2;
        dh *= 2;
        res2[i] = allocMem_dev(m, n, k);
        err2[i] = allocMem_dev(m, n, k);
        dimsinfo[i + 1] = {m, n, k, dh};
    }

    switch (DIM) {
            // case (1):
            //     for (i = 0; i < n_levels; i++) {
            //         innersinfo[i].istart = 0;     innersinfo[i].iend = 1;
            //         innersinfo[i].jstart = 0;     innersinfo[i].jend = 1;
            //         innersinfo[i].kstart = 1;     innersinfo[i].kend = dimsinfo[i].K - 1;
            //     }
            //     MGIterSolver = MGgaussseidel_1D;
            //     break;
            case (2):
                for (i = 0; i < n_levels; i++) {
                    innersinfo[i].istart = 0;     innersinfo[i].iend = 1;
                    innersinfo[i].jstart = 1;     innersinfo[i].jend = dimsinfo[i].N - 1;
                    innersinfo[i].kstart = 1;     innersinfo[i].kend = dimsinfo[i].K - 1;
                }
                // MGIterSolver = MGgaussseidel_2D;
                break;
            case (3):
                for (i = 0; i < n_levels; i++) {
                    innersinfo[i].istart = 1;     innersinfo[i].iend = dimsinfo[i].M - 1;
                    innersinfo[i].jstart = 1;     innersinfo[i].jend = dimsinfo[i].N - 1;
                    innersinfo[i].kstart = 1;     innersinfo[i].kend = dimsinfo[i].K - 1;
                }
                // MGIterSolver = MGgaussseidel_3D;
                break;
    }

    for (i = 0; i < n_levels; i++) {
        printf("%d_%d_%d_%d_%d_%d\n%d_%d_%d_%e\n", innersinfo[i].istart, innersinfo[i].iend,
                                                   innersinfo[i].jstart, innersinfo[i].jend,
                                                   innersinfo[i].kstart, innersinfo[i].kend,
                                                   dimsinfo[i].M, dimsinfo[i].N, dimsinfo[i].K, dimsinfo[i].DH);
    }

    dim3 block(32, 16, 1);
    dim3 grid;
    grid.x = (innersinfo[0].kend - innersinfo[0].kstart + block.x - 1) / block.x;
    grid.y = (innersinfo[0].jend - innersinfo[0].jstart + block.y - 1) / block.y;
    grid.z = (innersinfo[0].iend - innersinfo[0].istart + block.z - 1) / block.z;

    checkCudaErrors( cudaMemcpyToSymbol(dimsinfo_dev, dimsinfo, sizeof(DIMINFO)*n_levels) );
    checkCudaErrors( cudaMemcpyToSymbol(innersinfo_dev, innersinfo, sizeof(INNERS)*n_levels) );
    cublasHandle_t handle;
    checkCudaErrors( cublasCreate(&handle) );

    int iter = 0;
    double eps = 1.0e-6;
    double delta = 1.0;
    std::vector<double> delta_vector;
    type_t *delta_res;
    delta_res = allocMem_dev(M, N, K);

    double time_start, time_inter, time_end;
    time_start = seconds();
    // thrust::universal_vector<double> delta_vector;
    calcuDelta_dev <<<grid, block>>> (u, f, delta_res, 0);
    checkCudaErrors( cublasDnrm2(handle, M*N*K, delta_res, 1, &delta) );
    delta_vector.push_back(delta);
    time_inter = seconds();
    printf("iter: %6d\tdelta=%15g\ttime:%fs\n", iter, delta, time_inter-time_start);

    // int loop;
    // for (loop = 0; loop < 100; loop++) {
    while (delta > eps) {
        vcycle(u, f, 0);
        calcuDelta_dev <<<grid, block>>> (u, f, delta_res, 0);
        checkCudaErrors( cublasDnrm2(handle, M*N*K, delta_res, 1, &delta) );
        delta_vector.push_back(delta);

        // iter += n_cycles;
        iter += (n_cycles) * (n_levels > 1 ? (2 * n_levels - 1) : 1);
        time_inter = seconds();
        printf("iter: %6d\tdelta=%15g\ttime:%fs\n", iter, delta, time_inter-time_start);
        fflush(stdout);
    }

    time_end = seconds();
    printf("\nElapsed time: %fs\n", time_end-time_start);
    writeResult(u, delta_vector, M, N, K, DIM, time_end-time_start);

    checkCudaErrors( cudaFree(delta_res) );
    for (i = 0; i < n_levels - 1; i++) {
        // Only for 2D
        checkCudaErrors( cudaFree(res[i] ) );
        checkCudaErrors( cudaFree(res2[i]) );
        checkCudaErrors( cudaFree(err2[i]) );
    }
    
    return;

}



int main(int argc, char **argv)
{
    int M, N, K;
    int DIM;
    int method;
    // K: fast axis, N, M
    paramsCheck(argc, argv, &M, &N, &K, &DIM, &method, iterMth);

    size_t num = M * N * K;
    printf("num=%d\n", num);
    printf("sizeof=%d\n", sizeof(type_t));

    type_t *u_host = (type_t *)malloc(sizeof(type_t) * num);
    type_t *f_host = (type_t *)malloc(sizeof(type_t) * num);
    // type_t *r = (type_t *)malloc(sizeof(type_t) * num);
    // memset(r, 0, sizeof(type_t) * num);

    double len = 1.0;
    double DH = len / (K - 1);

    initialize(u_host, f_host, M, N, K, DIM, DH);
    // printf("%s, %d\n", __FILE__, __LINE__);

    type_t *u_dev, *f_dev;
    checkCudaErrors( cudaMalloc(&u_dev, sizeof(type_t) * num) );
    checkCudaErrors( cudaMalloc(&f_dev, sizeof(type_t) * num) );
    checkCudaErrors( cudaMemcpy(u_dev, u_host, sizeof(type_t) * num, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(f_dev, f_host, sizeof(type_t) * num, cudaMemcpyHostToDevice) );

    // multigrid method
    multigrid(u_dev, f_dev, M, N, K, DIM, DH);
    
    checkCudaErrors( cudaFree(u_dev) );
    checkCudaErrors( cudaFree(f_dev) );
    free(u_host);
    free(f_host);

    return 0;

}