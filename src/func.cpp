#include "libinclude.h"
#include "const.h"
#include "fundec.h"

// Function to calculate the dot product of two vectors of size n
template<typename T>
double dotProduct(T v1, T v2, size_t n){
    double result = 0;
    for (size_t i = 0; i < n; ++i) {
        result += v1[i] * v2[i];
    }
    return result;
}
template double dotProduct<float*>(float* v1, float* v2, size_t size);
template double dotProduct<double*>(double* v1, double* v2, size_t size);

// Function to calculate the cross product of two vectors and store it in the "out" vector
template<typename T1>
void crossProduct(T1 v_A, T1 v_B, double *out){
   out[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   out[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   out[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}
template void crossProduct<float*>(float*, float*, double*);

double F_0(double val){
    double a = 1-exp(-val*val);
    if(val>3.6){
        double e = exp(2*val*sqrt(M_PI)*log(2));
        double b = (e-1)/(e+1);
        return a - sqrt(M_PI)*val*b;
    }
    double b = sqrt(M_PI)*val*erf(val);
    return a-b;
}

double F_0_New(double DelZ, double DelZ2, double beta){
    double a = (1-exp(-DelZ2*beta*beta))/beta;
    double b = sqrt(M_PI)*DelZ*erf(beta*DelZ);
    return a-b;
}

double F_kl(double *ri, double *rj, double sigma, double psi, double beta, bool same_r, float **box){
    if(!same_r){
        double delX=ri[0]-rj[0];
        double delY=ri[1]-rj[1];
        double delZ=ri[2]-rj[2];
        sigma*=2*M_PI;
        psi*=2*M_PI;
        double norm_sigma_psi=sqrt(sigma*sigma+psi*psi);
        double a = cos(sigma*delX+psi*delY)/norm_sigma_psi;
        double b = exp(delZ*norm_sigma_psi)*erfc(beta*delZ+(norm_sigma_psi/(2*beta)))+exp(-delZ*norm_sigma_psi)*erfc(-beta*delZ+(norm_sigma_psi/(2*beta)));
        return a*b;
    }
    else{
        double norm_sigma_psi=2*M_PI*sqrt(sigma*sigma+psi*psi);
        return (2/norm_sigma_psi)*erfc(norm_sigma_psi/(2*beta));
    }
    return 0;
}

double M_n(double u, int n){
    if(n<2)return 0;
    else if(n==2){
        if(u<0 || u>n) return 0;
        else{
            return 1-abs(u-1);
        }
    }
    else{
        if(u<0 || u>n) return 0;
        else{
            return (u*M_n(u,n-1)/(n-1))+((n-u)*M_n(u-1,n-1)/(n-1));
        }
    }
}

complex<double>Coeff(double v, double w){
    const complex<double> t(0, 1);
    complex<double> bi_mi=exp(t*v*(w-1));
    complex<double> denox;
    for (double f = 0; f < w-1; f++){
        denox+=M_n(f+1,w)*exp(t*f*v);
    }
    if (denox == complex<double>(0, 0)) return 0;
    bi_mi/=denox;
    return bi_mi;
}

template<typename T2>
double error(T2 a, T2 b){
    return abs(a-b);
}
template double error<double>(double, double);

template<typename T3>
double percentReduction(T3 newValue, T3 oldValue){
    return 100*abs(newValue-oldValue)/oldValue;
}
template double percentReduction<double>(double, double);

