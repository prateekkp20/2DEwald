#include "libinclude.h"
#include "const.h"
#include "fundec.h"

double real(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, double cutoff){
    double real_energy=0;
    omp_set_num_threads(thread::hardware_concurrency());
    #pragma omp parallel for simd schedule(runtime) reduction(+: real_energy)
        for (int i = 0; i < natoms; i++){
            #pragma omp SIMD
            for (int j = 0; j < i; j++){
                if(i!=j){
                    double modR=dist(PosIons,i,j,box);
                    if(modR>cutoff)continue;
                    real_energy+=(ion_charges[i]*ion_charges[j]*erfc(betaa*modR))/modR;
                }
            }
        }
    
    return real_energy;
}