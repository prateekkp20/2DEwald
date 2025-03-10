#include "libinclude.h"
#include "const.h"
#include "fundec.h"
#include "header.h"

#define ENABLE_OMP 1
const complex<double> t(0,1);

double reciprocal_n2(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, int *K){
    double reciprocal_energy_i=0;
    #if defined ENABLE_OMP
        omp_set_num_threads(thread::hardware_concurrency());
    #endif
    // this is the loop for Ui
    double Length[3]={box[0][0],box[1][1],box[2][2]};
    for (int k = -K[0]; k < K[0]+1; k++){
        for (int l = -K[1]; l < K[1]+1; l++){
            if((k==0) && (l==0)) continue;
            double sigma= 2*M_PI*k*G[0][0];
            double psi= 2*M_PI*l*G[1][1];
            double norm_sigma_psi=sqrt(sigma*sigma+psi*psi);
            // The double loop sum over natoms in divided into two for loops:
            // (1) j<i
            // (2) j==i
            #if defined ENABLE_OMP
                #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_i)
            #endif
            // (1)
            for (int i = 0; i < natoms; i++){
                for (int j = 0; j < i; j++){
                    double delX=PosIons[3*i]-PosIons[3*j];
                    double delY=PosIons[3*i+1]-PosIons[3*j+1];
                    double delZ=PosIons[3*i+2]-PosIons[3*j+2];
                    reciprocal_energy_i+=2*ion_charges[i]*ion_charges[j]*F_kl_I(delX, delY, delZ, sigma, psi, norm_sigma_psi,betaa,box);
                }
            }

            #if defined ENABLE_OMP
                #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_i)
            #endif
            // (2)
            for (int i = 0; i < natoms; i++){
                reciprocal_energy_i+=ion_charges[i]*ion_charges[i]*F_kl_0(norm_sigma_psi,betaa);
            }
        }
    }
    
    return (M_PI/(2*Length[0]*Length[1]))*reciprocal_energy_i;
}

/*This Method is for the correction methods that we have introduced without the use of any particle-mesh methods*/
double reciprocal_modified(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, int *K){
    double reci_energy=0;
    #if defined ENABLE_OMP
        omp_set_num_threads(thread::hardware_concurrency());
        #pragma omp parallel for schedule(runtime) reduction(+: reci_energy) collapse(3)
    #endif
    for (int kx = -K[0]; kx < K[0]+1; kx++){
        for (int ky = -K[1]; ky < K[1]+1; ky++){
            for (int kz = -K[2]; kz < K[2]+1; kz++){
                if((kx==0) && (ky==0) && (kz==0))continue;

                int ii,jj,kk;
                if(kx<0) ii=(2*K[0]+1)+kx;
                else ii=kx;
                if(ky<0) jj=(2*K[1]+1)+ky;
                else  jj=ky;
                if(kz<0) kk=(2*K[2]+1)+kz;
                else  kk=kz;
                int temp=ii * ((2*K[2]+1) * (2*K[1]+1)) + jj * (2*K[2]+1) + kk;

                complex<double> sg=0;
                for (int  i = 0; i < natoms; i++){
                    double G_dot_r=2*M_PI*(kx*G[0][0]*PosIons[3*i]+ky*G[1][1]*PosIons[3*i+1]+kz*G[2][2]*PosIons[3*i+2]);
                    sg+=ion_charges[i]*(cos(G_dot_r)+t*sin(G_dot_r));
                }
                double norm_sg = norm(sg);
                
                //update energy
                reci_energy+=ExpFactor[temp]*norm_sg;
            }
        }
    }
    double L1 = box[0][0];
    double L2 = box[1][1];
    double L3 = box[2][2];
    return reci_energy*sqrt(M_PI)/(L1*L2);
}