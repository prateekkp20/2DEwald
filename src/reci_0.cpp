#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define ENABLE_OMP 1

double reci0(double **PosIons, float *ion_charges, int natoms, double betaa, float **box){
    double energy = 0;
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};

    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime) reduction(+: energy)
    #endif
    
    for (int  i = 0; i < natoms; i++){
        for (int j = 0; j < i; j++){
            energy+=ion_charges[i]*ion_charges[j]*F_0((PosIons[i][2]-PosIons[j][2])*betaa);
        }
    }
    
    return 2*sqrt(M_PI)*energy/(betaa*Length[0]*Length[1]);;
}