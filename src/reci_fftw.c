#include "libinclude.h"
#include "const.h"
#include "fundec.h"

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

double reciprocal_fftw(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K, int Grid, int n){
    // n: order of b-spline interpolation
    // initializing the new variables
    double G[3][3], m[3];
    double **u,**x_direc, **y_direc;
    u= new double * [natoms];
    x_direc = new double * [natoms];
    y_direc = new double * [natoms];

    for (int  i = 0; i < natoms; i++){
        u[i] = new double  [2];
        x_direc[i] = new double  [K];
        y_direc[i] = new double  [K];
    }

    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    int n_max=2;

    omp_set_num_threads(thread::hardware_concurrency());
    // this is for Ui
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

    gsl_function F;
    // F.function = &integrand_reciprocal; // Set the function to integrate
    // reciprocal_n_params params = {PosIons, ion_charges, natoms, betaa, box, K};
    // F.params = &params;

    double result, error;
    // gsl_integration_qagi(&F, 1e-7, 1e-2, 1000, workspace, &result, &error);
    // gsl_integration_workspace_free(workspace); // Free workspace memory
    
    double reciprocal_energy_o=0;

    // this is the loop for Uo
    #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_o) collapse(2)
    for (int  i = 0; i < natoms; i++){
        for (int j = 0; j < natoms; j++){
            reciprocal_energy_o+=ion_charges[i]*ion_charges[j]*F_0(PosIons[i][2]-PosIons[j][2],betaa);
        }
    }
    return sqrt(M_PI)*reciprocal_energy_o/(Length[0]*Length[1])+result;
}