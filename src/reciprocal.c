#include "libinclude.h"
#include "const.h"
double F_kl(double *r, double i,double beta, ){
    double F_kl;

    return F_kl;
}

double reciprocal_n2(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K){
    double reciprocal_energy=0;
    for (int k = -K; k < K+1; k++){
        for (int l = -K; l < K+1; l++){
            if((k==0) && (l==0)) continue;
            for (int i = 0; i < natoms; i++){
                for (int j = 0; j < natoms; j++){
                    reciprocal_energy+=ion_charges[i]*ion_charges[j]*F_kl(PosIons[i]-PosIons[j])
                }
                
            }
            
        }
    }
    return reciprocal_energy;
}