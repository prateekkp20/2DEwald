#include "libinclude.h"
#include "const.h"
#include "fundec.h"

double reciprocal_n2(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K){
    double reciprocal_energy_i=0;
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
                    reciprocal_energy_i+=2*ion_charges[i]*ion_charges[j]*F_kl(PosIons[i],PosIons[j],sigma,psi,betaa,false,box);
                }
            }
            #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_i)
            for (int i = 0; i < natoms; i++){
                reciprocal_energy_i+=ion_charges[i]*ion_charges[i]*F_kl(PosIons[i],PosIons[i],sigma,psi,betaa,true,box);
            }
        }
    }
    
    return (M_PI/(2*Length[0]*Length[1]))*reciprocal_energy_i;
}