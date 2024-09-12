#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define ENABLE_OMP 1

double real(double *PosIons2, double *charge_prod, int natoms, double betaa, double **box, double cutoff){
    double real_energy=0;

    #if defined ENABLE_OMP
        omp_set_num_threads(thread::hardware_concurrency());
        #pragma omp parallel for simd schedule(runtime) reduction(+: real_energy)
    #endif

        for (int i = 1; i < natoms; i++){
            for (int j = 0; j < i; j++){
                    double modR=dist(PosIons2,i,j,box);
                    if(modR>cutoff)continue;
                    real_energy+=(charge_prod[i*(i-1)/2+j]*erfc(betaa*modR))/modR;
            }
        }
    
    return real_energy;
}