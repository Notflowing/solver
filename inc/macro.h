#ifndef __MACRO__
#define __MACRO__

typedef struct PARAMS {
    int M;
    int N;
    int K;
    int DIM;
    int method;
    double DH;
} PARAMS;

typedef struct INNERS {
    int istart;
    int iend;
    int jstart;
    int jend;
    int kstart;
    int kend;
} INNERS;

#define RELAX 1.25

inline double seconds()
{
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

inline void swapPtr(type_t **a, type_t **b)
{
    type_t *temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

#endif