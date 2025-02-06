#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define ENABLE_OMP 1

double reciprocal_n2(double **PosIons, double *charge_prod, double *ion_charges, int natoms, double betaa, double **box, int *K){
    double reciprocal_energy_i=0;
    #if defined ENABLE_OMP
        omp_set_num_threads(thread::hardware_concurrency());
    #endif
    // this is the loop for Ui
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    for (int k = -K[0]; k < K[0]+1; k++){
        for (int l = -K[1]; l < K[1]+1; l++){
            if((k==0) && (l==0)) continue;
            double sigma= k/Length[0];
            double psi= l/Length[1];
            // The double loop sum over natoms in divided into two for loops:
            // (1) j<i
            // (2) j==i
            #if defined ENABLE_OMP
                #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_i)
            #endif
            // (1)
            for (int i = 0; i < natoms; i++){
                for (int j = 0; j < i; j++){
                    reciprocal_energy_i+=2*charge_prod[i*(i-1)/2+j]*F_kl_I(PosIons[i],PosIons[j],sigma,psi,betaa,box);
                }
            }

            #if defined ENABLE_OMP
                #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_i)
            #endif
            // (2)
            for (int i = 0; i < natoms; i++){
                reciprocal_energy_i+=ion_charges[i]*ion_charges[i]*F_kl_0(sigma,psi,betaa);
            }
        }
    }
    
    return (M_PI/(2*Length[0]*Length[1]))*reciprocal_energy_i;
}