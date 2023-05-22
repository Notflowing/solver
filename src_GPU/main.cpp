#include "header.h"

// std::map<std::string, int> iterMth = {{1, "JACOBI"}, {2, "GAUSSSEIDEL"}, {3, "SOR"}, {4, "RBGS"}, {5, "MGRBGS"}};
const char *iterMth[] = {"DEFAULT", "JACOBI", "GAUSSSEIDEL", "SOR", "RBGS", "MGRBGS"};
__constant__ int dim_dev[3];
__constant__ int inn_dev[6];
cublasHandle_t handle;


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
    type_t *u_cur = (type_t *)malloc(sizeof(type_t) * num);
    type_t *u_new = (type_t *)malloc(sizeof(type_t) * num);
    type_t *error = (type_t *)malloc(sizeof(type_t) * num);
    type_t *f = (type_t *)malloc(sizeof(type_t) * num);

    double len = 1.0;
    double DH = len / (K - 1);
    initU(u_cur, u_new, error, M, N, K);
    initF(f, M, N, K, DIM, DH);

    double eps = 1.0e-4;
    double delta = 1.0;
    int iter = 0;
    std::vector<double> delta_vector;

    PARAMS params = {M, N, K, DIM, method, DH};
    INNERS inners;
    int methodRank;
    configstcl(&inners, &methodRank, &params);

    // GPU configuration
    type_t *u_cur_dev;
    type_t *u_new_dev;
    type_t *error_dev;
    type_t *f_dev;
    checkCudaErrors( cudaMalloc((void **)&u_cur_dev, sizeof(type_t) * num) );
    checkCudaErrors( cudaMalloc((void **)&u_new_dev, sizeof(type_t) * num) );
    checkCudaErrors( cudaMalloc((void **)&error_dev, sizeof(type_t) * num) );
    checkCudaErrors( cudaMalloc((void **)&f_dev, sizeof(type_t) * num) );
    checkCudaErrors( cudaMemcpy(u_cur_dev, u_cur, sizeof(type_t) * num, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(u_new_dev, u_new, sizeof(type_t) * num, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(error_dev, error, sizeof(type_t) * num, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(f_dev, f, sizeof(type_t) * num, cudaMemcpyHostToDevice) );
    type_t *u_cur_host = u_cur;
    free(u_new);
    free(error);
    free(f);
    u_cur = u_cur_dev;
    u_new = u_new_dev;
    error = error_dev;
    f = f_dev;
    int dim[3] = {M, N, K};
    int inn[6] = {inners.istart, inners.iend, inners.jstart, inners.jend, inners.kstart, inners.kend};
    checkCudaErrors( cudaMemcpyToSymbol(dim_dev, dim, sizeof(int)*3) );
    checkCudaErrors( cudaMemcpyToSymbol(inn_dev, inn, sizeof(int)*6) );
    checkCudaErrors( cublasCreate(&handle) );


    double time_start, time_inter, time_end;
    time_start = seconds();
    // int loop;
    // for (loop = 0; loop < 20000; loop++) {
    while (delta > eps) {
        // iterative method
        delta = iterMethod(u_cur, u_new, error, f, params, inners, methodRank);
        delta_vector.push_back(delta);
        swapPtr(&u_cur, &u_new);
        if (iter % 1000 == 0) {
            time_inter = seconds();
            printf("iter: %6d\tdelta=%g\ttime:%fs\n", iter, delta, time_inter-time_start);
            fflush(stdout);
        }
        iter++;
    }

    checkCudaErrors( cublasDestroy(handle) );
    checkCudaErrors( cudaMemcpy(u_cur_host, u_cur, sizeof(type_t) * num, cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaFree(u_cur) );
    checkCudaErrors( cudaFree(u_new) );
    checkCudaErrors( cudaFree(error) );
    checkCudaErrors( cudaFree(f) );
    u_cur = u_cur_host;

    writeResult(u_cur, delta_vector, params);
    free(u_cur);
    time_end = seconds();
    printf("Elapsed time: %fs\n", time_end-time_start);
    
    return 0;
}

