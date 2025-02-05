// In this file I would be doing the testing for the 2D Ewald with a hat function correction, and the various functions needed to calculate the correction factor 
#include "libinclude.h"
#include "const.h"
#include "fundec.h"
#include "header.h"

struct C_params{
    double t; double L; double gamm; double n;
};

struct constant_params{
    int kx; int ky; int kz; double lx; double ly; double lz; double gamma;
};

double hat_function(double s, double gamm, double L){
    double result = 1/(1+exp(-gamm*(0.5*L+s))) + 1/(1+exp(-gamm*(0.5*L-s))) - 1;
    return result;
}

double C_integrand(double s, void *params){
    C_params* p = (struct C_params*)params;
    double t = p->t;
    double L = p->L;
    double gamm = p->gamm;
    double n = p->n;

    double result=exp(-s*t*s*t)*cos(2*M_PI*n*s/L)*hat_function(s,gamm,L);
    return result;
}

double C(double t, double n, double gamm, double L){

    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(200);
    gsl_function F;
    F.function = &C_integrand; // Set the function to integrate
    C_params params = {t, L, gamm, n};
    F.params = &params;
    double result, error;
    gsl_integration_qagi(&F, 1e-4, 1e-2, 200, workspace, &result, &error);
    gsl_integration_workspace_free(workspace); // Free workspace memory

    return result/L;
}

double constant_integrand(double t, void *params){
    constant_params* p = (struct constant_params*)params;
    int kx = p->kx; int ky = p->ky; int kz = p->kz; double lx = p->lx; double ly = p->ly; double lz = p->lz; double gamma = p->gamma;
    double g = M_PI*M_PI*((kx*kx)/(lx*lx)+(ky*ky)/(ly*ly));
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
    // gsl_integration_qng(&F, 0, beta, 1e-4, 1e-2, &result, &error, &neval);
    // int gsl_integration_qng(const gsl_function *f, double a, double b, double epsabs, double epsrel, double *result, double *abserr, size_t *neval)

    gsl_integration_qag(&F,0, beta, 1e-4, 1e-2, 200, 3, workspace, &result, &error);// this has less error 
    // int gsl_integration_qag(const gsl_function *f, double a, double b, double epsabs, double epsrel, size_t limit, int key, gsl_integration_workspace *workspace, double *result, double *abserr)

    cout<<"Error: "<<error<<endl;
    gsl_integration_workspace_free(workspace); // Free workspace memory

    return result;
}

int main(){
    cout<<fixed<<setprecision(15)<<C(20,3,400,1)<<endl;
    cout<<fixed<<setprecision(15)<<constantterm(2,2,2,10,10,10,5.42/10,400)<<endl;
    cout<<"FILE RUNS"<<endl;
    return 0;
}