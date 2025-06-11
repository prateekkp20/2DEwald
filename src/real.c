#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define ENABLE_OMP 1

double real(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, double cutoff){
    double real_energy=0;

    #if defined ENABLE_OMP
        #pragma omp parallel for simd schedule(runtime) reduction(+: real_energy)
    #endif

        for (int i = 0; i < natoms; i++){
            for (int j = 0; j < i; j++){
                    double modR=dist(PosIons,i,j,box);
                    if(modR>cutoff)continue;
                    real_energy+=(ion_charges[i]*ion_charges[j]*erfc(betaa*modR))/modR;
            }
        }
    
    return real_energy;
}