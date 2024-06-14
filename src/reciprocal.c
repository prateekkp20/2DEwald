#include "libinclude.h"
#include "const.h"
#include "fundec.h"

double F_0(double Del_z, double beta){
    double a = (1-exp(-Del_z*Del_z*beta*beta))/beta;
    double b = sqrt(M_PI)*Del_z*erf(beta*Del_z);
    return a-b;
}

double F_kl(double *ri, double *rj, double sigma, double psi, double beta, bool same_r){
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

double reciprocal_n2(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K){
    double reciprocal_energy_i=0;
    double reciprocal_energy_o=0;
    omp_set_num_threads(thread::hardware_concurrency());
    // this is the loop for Ui
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    for (int k = -K; k < K+1; k++){
        for (int l = -K; l < K+1; l++){
            if((k==0) && (l==0)) continue;
            double sigma= k/Length[0];
            double psi= l/Length[1];
            // The double loop sum over natoms in divided into two for loops:
            // (1) j<i
            // (2) j==i
            #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_i)
            for (int i = 0; i < natoms; i++){
                for (int j = 0; j < i; j++){
                    reciprocal_energy_i+=2*ion_charges[i]*ion_charges[j]*F_kl(PosIons[i],PosIons[j],sigma,psi,betaa,false);
                }
            }
            #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_i)
            for (int i = 0; i < natoms; i++){
                reciprocal_energy_i+=ion_charges[i]*ion_charges[i]*F_kl(PosIons[i],PosIons[i],sigma,psi,betaa,true);
            }
        }
    }
    #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_o) collapse(2)
    // this is the loop for Uo
    for (int  i = 0; i < natoms; i++){
        for (int j = 0; j < natoms; j++){
            reciprocal_energy_o+=ion_charges[i]*ion_charges[j]*F_0(PosIons[i][2]-PosIons[j][2],betaa);
        }
    }
    
    return sqrt(M_PI)*reciprocal_energy_o/(Length[0]*Length[1]) +(M_PI/(2*Length[0]*Length[1]))*reciprocal_energy_i;
}