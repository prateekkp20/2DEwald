#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define ENABLE_OMP 1

vector<double> realnreci0(double *PosIons2, double *charge_prod, int natoms, double betaa, double **box, double cutoff){
    // This term will contain the real energy part as energy[0] and reciprocal energy k=0 as energy[1]
    vector<double> energy(2,0);
    double energy0=0,energy1=0;
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    #if defined ENABLE_OMP
        omp_set_num_threads(thread::hardware_concurrency());
        #pragma omp parallel for simd schedule(runtime) reduction(+: energy0,energy1)
    #endif
    for (int i = 1; i < natoms; i++){
        for (int j = 0; j < i; j++){

            // real energy
            double modR=dist(PosIons2,i,j,box);
            if(modR<=cutoff)
                energy0+=(charge_prod[i*(i-1)/2+j]*erfc(betaa*modR))/modR;
            
            // reciprocal energy k=0
            energy1+=charge_prod[i*(i-1)/2+j]*F_0((PosIons2[i*3+2]-PosIons2[j*3+2])*betaa);
        }
    }
    energy1*=2*sqrt(M_PI)/(betaa*Length[0]*Length[1]);
    energy[0]=energy0;
    energy[1]=energy1;
    return energy;
}