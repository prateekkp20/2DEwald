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
template void crossProduct<double*>(double*, double*, double*);

double F_0(double val){
    val=fabs(val); //making the input positive
    // this is the rational approximation of the error function, mentioned in Abramowitz and Stegun, Handbook of Mathematical functions
    // this give the max errror of e-7 wrt to the actual erf()
    // float a1 = 0.254829592,a2 = -0.284496736,a3 = 1.421413741,a4 = -1.453152027,a5 = 1.061405429; 
    double exp_x2 = exp(-val*val);
    double t, t1 =  t  = 1/(1+0.3275911*val);
    double erfx = 1 - exp_x2*(0.254829592*t - 0.284496736*(t*=t1) + 1.421413741*(t*=t1) - 1.453152027*(t*=t1) + 1.061405429*(t*=t1));
    double b = 1.772453850905515*val*erfx;
    return 1-exp_x2-b;
}

double F_kl(double *ri, double *rj, double sigma, double psi, double beta, bool same_r, double **box){
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

double F_kl_0(double sigma, double psi, double beta){
    double norm_sigma_psi=2*M_PI*sqrt(sigma*sigma+psi*psi);
    return (2/norm_sigma_psi)*erfc(norm_sigma_psi/(2*beta));
}

double F_kl_I(double *ri, double *rj, double sigma, double psi, double beta, double **box){
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

template<typename T4>
T4* linspace(T4 a, T4 b, T4 num){
    T4 * line  = new T4[b-a+1];
    for (int i = 0; i < b-a+1; i++){
        line[i] = a+i;
    }
    return line;
}
template int* linspace(int,int,int);