#ifndef __FUNCTIONS__
#define __FUNCTIONS__

void paramsCheck(int argc, char **argv, int *M, int *N, int *K, int *DIM,
                 int *method, const char **iterMth);
void initU(type_t *u_cur, type_t *u_new, const int Mdim, const int Ndim, const int Kdim);
void initF(type_t *f, const int Mdim, const int Ndim, const int Kdim);
void configstcl(INNERS *inners, iterStencil_t *iterstcl, PARAMS *params);
double iterMethod(type_t *u_cur, type_t *u_new, type_t *f, PARAMS params,
                  INNERS inners, iterStencil_t iterstcl);

void writeResult(type_t *u_cur, std::vector<double> &delta, PARAMS params);

// iteration method
void jacobi_1D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
               const double DH2, const int Mdim, const int Ndim, const int Kdim,
               const int i, const int j, const int k);
void jacobi_2D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
               const double DH2, const int Mdim, const int Ndim, const int Kdim,
               const int i, const int j, const int k);
void jacobi_3D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
               const double DH2, const int Mdim, const int Ndim, const int Kdim,
               const int i, const int j, const int k);

void gaussseidel_1D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
                    const double DH2, const int Mdim, const int Ndim, const int Kdim,
                    const int i, const int j, const int k);
void gaussseidel_2D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
                    const double DH2, const int Mdim, const int Ndim, const int Kdim,
                    const int i, const int j, const int k);
void gaussseidel_3D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
               const double DH2, const int Mdim, const int Ndim, const int Kdim,
               const int i, const int j, const int k);

void sor_1D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
            const double DH2, const int Mdim, const int Ndim, const int Kdim,
            const int i, const int j, const int k);
void sor_2D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
            const double DH2, const int Mdim, const int Ndim, const int Kdim,
            const int i, const int j, const int k);
void sor_3D(type_t * __restrict__ u_cur, type_t * __restrict__ u_new, type_t * __restrict__ f,
            const double DH2, const int Mdim, const int Ndim, const int Kdim,
            const int i, const int j, const int k);



#endif