#include "libinclude.h"
#include "const.h"
#include "fundec.h"

struct reciprocal_n_params {
  double** PosIons;
  float* ion_charges;
  int natoms;
  double betaa;
  float** box;
  int K;
};

double integrand_reciprocal(double h, void *params){
    // this integral is over the variable h;
    reciprocal_n_params* p = (struct reciprocal_n_params*)params;
    double **PosIons = p->PosIons;
    float *ion_charges = p->ion_charges;
    int natoms = p->natoms;
    double betaa = p->betaa;
    float **box = p->box;
    int K = p->K;

    double reciprocal_energy_i=0;
    //Length of the sides of the unit cell
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_i) collapse(2)
    for (int k = -K; k < K+1; k++){
        for (int l = -K; l < K+1; l++){
            if((k==0) && (l==0))continue;
            complex<double> s_kh=0;
            complex<double> t(0,1);
            double G[3]={2*M_PI*k/Length[0], 2*M_PI*l/Length[1], h};
            for (int  i = 0; i < natoms; i++){
                double G_dot_r=G[0]*PosIons[i][0]+G[1]*PosIons[i][1]+G[2]*PosIons[i][2];
                complex<double> charge(ion_charges[i],0.0);
                s_kh+=charge*(cos(G_dot_r)+t*sin(G_dot_r));
            }
            double norm_sg = norm(s_kh);
            double norm_g = G[0]*G[0]+G[1]*G[1]+G[2]*G[2];
            reciprocal_energy_i+=norm_sg/(norm_g*exp(norm_g/(4*betaa*betaa)));
        }
    }
    reciprocal_energy_i*=1/(Length[0]*Length[1]);
    return reciprocal_energy_i;
}

double reciprocal_kawata(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K) {
    omp_set_num_threads(thread::hardware_concurrency());
    // this is for Ui
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

    gsl_function F;
    F.function = &integrand_reciprocal; // Set the function to integrate
    reciprocal_n_params params = {PosIons, ion_charges, natoms, betaa, box, K};
    F.params = &params;

    double result, error;
    gsl_integration_qagi(&F, 1e-7, 1e-2, 1000, workspace, &result, &error);
    gsl_integration_workspace_free(workspace); // Free workspace memory

    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    double reciprocal_energy_o=0;

    // this is the loop for Uo
    #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_o) collapse(2)
    for (int  i = 0; i < natoms; i++){
        for (int j = 0; j < natoms; j++){
            reciprocal_energy_o+=ion_charges[i]*ion_charges[j]*F_0(PosIons[i][2]-PosIons[j][2],betaa);
        }
    }
    return sqrt(M_PI)*reciprocal_energy_o/(Length[0]*Length[1])+result;
    // return result;
}
