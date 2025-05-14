// This file has the functions, and data structures to calculate the Screening Function for our modified method 
#include "libinclude.h"
#include "const.h"
#include "fundec.h"
#include "header.h"

struct C_params{
    double t; double L; double gamma; double n;
};

struct constant_params{
    int kx; int ky; int kz; double lx; double ly; double lz; double gamma;
};

double tophat(double s, double gamma, double L){
    double result = 1/(1+exp(-gamma*(0.5*L+s))) + 1/(1+exp(-gamma*(0.5*L-s))) - 1; // sigmoid type function
    // double result = tanh(gamma*(s+0.5*L)) - tanh(gamma*(s-0.5*L));
    return result;
}

double C_integrand(double s, void *params){
    C_params* p = (struct C_params*)params;
    double t = p->t;
    double L = p->L;
    double gamma = p->gamma;
    double n = p->n;

    double result=exp(-s*t*s*t)*cos(2*M_PI*n*s/L)*tophat(s,gamma,L);
    return result;
}

double C(double t, double n, double gamma, double L){

    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(200);
    gsl_function F;
    F.function = &C_integrand; // Set the function to integrate
    C_params params = {t, L, gamma, n};
    F.params = &params;
    double result, error;
    gsl_integration_qagi(&F, 1e-8, 1e-6, 200, workspace, &result, &error);
    gsl_integration_workspace_free(workspace); // Free workspace memory

    return result/L;
}

double constant_integrand(double t, void *params){
    constant_params* p = (struct constant_params*)params;
    int kx = p->kx; int ky = p->ky; int kz = p->kz; double lx = p->lx; double ly = p->ly; double lz = p->lz; double gamma = p->gamma;
    double m[3];
    for (int t = 0; t < 3; t++){
        m[t]=kx*G[0][t]+ky*G[1][t];   
    }
    double g=M_PI*M_PI*dotProduct(m,m,3);
    double result = (1/(t*t))*C(t,kz,gamma,lz)*exp(-g/(t*t));
    return result;
}

double constantterm(int kx, int ky, int kz, double lx, double ly, double lz, double beta, double gamma){

    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(200);
    gsl_function F;
    F.function = &constant_integrand; // Set the function to integrate
    constant_params params = {kx, ky, kz, lx, ly, lz, gamma};
    F.params = &params;
    double result, error;
    size_t neval;
    gsl_integration_qng(&F, 0, beta, 1e-6, 1e-6, &result, &error, &neval);
    // int gsl_integration_qng(const gsl_function *f, double a, double b, double epsabs, double epsrel, double *result, double *abserr, size_t *neval)

    // gsl_integration_qag(&F,0, beta, 1e-4, 1e-2, 200, 3, workspace, &result, &error);// this has less error 
    // int gsl_integration_qag(const gsl_function *f, double a, double b, double epsabs, double epsrel, size_t limit, int key, gsl_integration_workspace *workspace, double *result, double *abserr)

    gsl_integration_workspace_free(workspace); // Free workspace memory

    return result;
}