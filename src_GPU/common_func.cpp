#include "header.h"

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



void writeResult(type_t *u_cur, std::vector<double> &delta, PARAMS params)
{
    int Mdim = params.M;
    int Ndim = params.N;
    int Kdim = params.K;
    int DIM = params.DIM;
    int method = params.method;
    double DH = params.DH;
    extern const char *iterMth[];

    // int i, j, k;
    // for (i = 0; i < Mdim; i++) {
    //     for (j = 0; j < Ndim; j++) {
    //         for (k = 0; k < Kdim; k++) {
    //             printf("%f\n", u_cur[INDEX(i, j, k)]);
    //         }
    //     }
    // }

    FILE *fp;
    char filename[128];
    sprintf(filename, "output/%s_%dD_GPU.bin", iterMth[method], DIM);
    fp = fopen(filename, "wb");
    fwrite(u_cur, sizeof(type_t), Mdim*Ndim*Kdim, fp);
    fclose(fp);

    FILE *fp1;
    char filename1[128];
    sprintf(filename1, "output/%s_%dD_dim_type_GPU.txt", iterMth[method], DIM);
    fp1 = fopen(filename1, "w");
    fprintf(fp1, "DIM %d %d %d sizeof %d", Mdim, Ndim, Kdim, sizeof(type_t));
    fclose(fp1);

    FILE *fp2;
    char filename2[128];
    sprintf(filename2, "output/%s_%dD_delta_GPU.bin", iterMth[method], DIM);
    fp2 = fopen(filename2, "wb");
    fwrite(delta.data(), sizeof(double), delta.size(), fp2);
    fclose(fp2);

}

void configstcl(INNERS *inners, int *rank, PARAMS *params)
{
    int Mdim = params->M;
    int Ndim = params->N;
    int Kdim = params->K;
    int DIM = params->DIM;
    int method = params->method;
    
    int istart, iend;
    int jstart, jend;
    int kstart, kend;

// #if defined(GPU)
//     thrust::universal_vector<iterStencil_t> ITER(9);
//     ITER[0] = jacobi_1D;
//     ITER[1] = jacobi_2D;
//     ITER[2] = jacobi_3D;
//     ITER[3] = gaussseidel_1D;
//     ITER[4] = gaussseidel_2D;
//     ITER[5] = gaussseidel_3D;
//     ITER[6] = sor_1D;
//     ITER[7] = sor_2D;
//     ITER[8] = sor_3D;
// #else

// #endif

    switch (DIM) {
        case (1):
            istart = 0;     iend = 1;
            jstart = 0;     jend = 1;
            kstart = 1;     kend = Kdim - 1;
            break;
        case (2):
            istart = 0;     iend = 1;
            jstart = 1;     jend = Ndim - 1;
            kstart = 1;     kend = Kdim - 1;
            break;
        case (3):
            istart = 1;     iend = Mdim - 1;
            jstart = 1;     jend = Ndim - 1;
            kstart = 1;     kend = Kdim - 1;
            break;
    }
    *inners = {istart, iend, jstart, jend, kstart, kend};

    *rank = (method-1) * 3 + DIM - 1;

    printf("%d\t%d\n%d\t%d\n%d\t%d\n%d\n", istart, iend, jstart, jend, kstart, kend, *rank + 1);
    fflush(stdout);

}