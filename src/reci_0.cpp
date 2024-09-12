#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define ENABLE_OMP 1

double reci0(double *PosIons2, double *charge_prod, int natoms, double betaa, double **box){
    double energy = 0;
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};

    #if defined ENABLE_OMP
        omp_set_num_threads(thread::hardware_concurrency());
        #pragma omp parallel for schedule(runtime) reduction(+: energy)
    #endif
    
    for (int  i = 1; i < natoms; i++){
        for (int j = 0; j < i; j++){
            energy+=charge_prod[i*(i-1)/2+j]*F_0((PosIons2[i*3+2]-PosIons2[j*3+2])*betaa);
        }
    }
    
    return 2*sqrt(M_PI)*energy/(betaa*Length[0]*Length[1]);;
}