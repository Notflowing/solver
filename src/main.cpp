#include "header.h"

// std::map<std::string, int> iterMth = {{1, "JACOBI"}, {2, "GAUSSSEIDEL"}, {3, "SOR"}, {4, "RBGS"}, {5, "MGRBGS"}};
const char *iterMth[] = {"DEFAULT", "JACOBI", "GAUSSSEIDEL", "SOR", "RBGS", "MGRBGS"};

int main(int argc, char **argv)
{
    int M, N, K;
    int DIM;
    int method;
    paramsCheck(argc, argv, &M, &N, &K, &DIM, &method, iterMth);

    size_t num = M * N * K;
    // printf("num=%d\n", num);
    // printf("sizeof=%d\n", sizeof(type_t));
    type_t *u_cur = (type_t *)malloc(sizeof(type_t) * num);
    type_t *u_new = (type_t *)malloc(sizeof(type_t) * num);
    type_t *f = (type_t *)malloc(sizeof(type_t) * num);

    double len = 1.0;
    double DH = len / (K - 1);
    initU(u_cur, u_new, M, N, K);
    initF(f, M, N, K);

#if defined(__CUDA_ARCH__)
    type_t *u_cur_dev;
    type_t *u_new_dev;
    type_t *f_dev;
    checkCudaErrors( cudaMalloc((void **)u_cur_dev, sizeof(type_t) * num) );
#endif

    double eps = 1.0e-6;
    double delta = 1.0;
    int iter = 0;

    PARAMS params = {M, N, K, DIM, method, DH};
    INNERS inners;
    iterStencil_t iterstcl;
    configstcl(&inners, &iterstcl, &params);

    std::vector<double> delta_vector;
    double time_start, time_inter, time_end;
    time_start = seconds();
    while (delta > eps) {
        // iterative method
        delta = iterMethod(u_cur, u_new, f, params, inners, iterstcl);
        delta_vector.push_back(delta);
        swapPtr(&u_cur, &u_new);
        if (iter % 1000 == 0) {
            time_inter = seconds();
            printf("iter: %6d\tdelta=%g\ttime:%fs\n", iter, delta, time_inter-time_start);
            fflush(stdout);
        }
        iter++;
    }

    writeResult(u_cur, delta_vector, params);
    time_end = seconds();
    printf("Elapsed time: %fs\n", time_end-time_start);
    
    return 0;
}

