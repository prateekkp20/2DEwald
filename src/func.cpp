#include "libinclude.h"
#include "const.h"
#include "fundec.h"

template<typename T>
double dotProduct(T v1, T v2, size_t size){
    double result = 0;
    for (size_t i = 0; i < size; ++i) {
        result += v1[i] * v2[i];
    }
    return result;
}

template double dotProduct<float*>(float* v1, float* v2, size_t size);
template double dotProduct<double*>(double* v1, double* v2, size_t size);

template<typename T1>
void crossProduct(T1 v_A, T1 v_B, double *out){
   out[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   out[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   out[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

template void crossProduct<float*>(float*, float*, double*);