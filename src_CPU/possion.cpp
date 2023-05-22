#include "header.h"

void init_u_host(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ error,
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
                error[INDEX(i, j, k)] = 0.0;

            }
        }
    }
}

void initU(type_t *u_cur, type_t *u_new, type_t *error, const int Mdim, const int Ndim, const int Kdim)
{
    init_u_host(u_cur, u_new, error, Mdim, Ndim, Kdim);
}



void init_f_host(type_t * __restrict__ f, const int Mdim, const int Ndim, const int Kdim, const int DIM, const double DH)
{
    int i, j, k;
    long long idx;
    // auto func = [](const int i, const int j, const int k, const double DH){return sin(KHZ * PI * k * DH);};
    std::function<type_t(const int, const int, const int, const double)> func;
    switch (DIM) {
        case (1):
            func = [](const int i, const int j, const int k, const double DH) \
                   {return sin(KHZ * PI * k * DH);};
            break;
        case (2):
            func = [](const int i, const int j, const int k, const double DH) \
                   {return sin(KHZ * PI * k * DH) * sin(KHZ * PI * j * DH);};
            break;
        case (3):
            func = [](const int i, const int j, const int k, const double DH) \
                   {return sin(KHZ * PI * k * DH) * sin(KHZ * PI * j * DH) * sin(KHZ * PI * i * DH);};
            break;
    }
    
    for (i = 0; i < Mdim; i++) {
        for (j = 0; j < Ndim; j++) {
            for (k = 0; k < Kdim; k++) {
                idx = INDEX(i, j, k);
                // f[idx] = 0.0;
                // if (i == Mdim / 2 && j == Ndim / 2 && k == Kdim / 2) {
                //     f[idx] = -(Kdim-1);
                // }
                // f[idx] = sin(KHZ * PI * k * DH) * sin(KHZ * PI * j * DH);
                f[idx] = - 2 * KHZ * KHZ * PI * PI * func(i, j, k, DH);

            }
        }
    }
    
}

void initF(type_t *f, const int Mdim, const int Ndim, const int Kdim, const int DIM, const double DH)
{

    init_f_host(f, Mdim, Ndim, Kdim, DIM, DH);

}