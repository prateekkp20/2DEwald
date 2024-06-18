#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define REAL 0
#define IMAG 1

double reciprocal_pppm(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K, int Grid[], int n[]){
    /* Side lengths of the unit cell */
    omp_set_num_threads(thread::hardware_concurrency());
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))}; 
    double reciprocal_energy_o=0;
    double reciprocal_energy_i=0.0;
    double G[3][3];// Reciprocal Vector
    double **u; // Fractional coordinates of the atoms in the unit cell
    double **x_direc, **y_direc, **z_direc;
    u= new double * [natoms];
    x_direc = new double * [natoms];
    y_direc = new double * [natoms];
    z_direc = new double * [natoms];
    for (int  i = 0; i < natoms; i++){
        u[i] = new double  [2]; // We only need these in x and y direction 
        x_direc[i] = new double  [Grid[0]];
        y_direc[i] = new double  [Grid[1]];
        z_direc[i] = new double  [Grid[2]];
    }

    // Volume Calculations
    double A[3];
    double C[3]={box[2][0],box[2][1],box[2][2]};
    crossProduct(box[0],box[1],A);
    double volume = dotProduct(A,C,3);

    // Calculating the reciprocal vectors
    crossProduct(box[1],box[2],G[0]);
    crossProduct(box[2],box[0],G[1]);
    crossProduct(box[0],box[1],G[2]);
    for (int x = 0; x < 3; x++)
        for (int q = 0; q < 3; q++)
            G[x][q] /= volume;

    // Calculating the fractional coordinates in x and y directions
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 2; j++){
            u[i][j]=Grid[j]*dotProduct(PosIons[i],G[j],3);
        }
    }
    #pragma omp parallel for simd
    // Calculating the cofficients in the x,y and z directions for the Q Matrix
    for (int i = 0; i < natoms; i++){
        // for X direction
        for (int  k1 = 0; k1 < Grid[0]; k1++){
            x_direc[i][k1]=0;
            for (int  n1 = ceil((u[i][0]-k1-n[0])/Grid[0]); n1 <= floor((u[i][0]-k1)/Grid[0]); n1++){
                x_direc[i][k1]+=M_n(u[i][0]-k1-n1*Grid[0],n[0]);
            }
        }
        // for Y direction
        for (int  k2 = 0; k2 < Grid[1]; k2++){
            y_direc[i][k2]=0;
            for (int  n2 = ceil((u[i][1]-k2-n[1])/Grid[1]); n2 <= floor((u[i][1]-k2)/Grid[1]); n2++){
                y_direc[i][k2]+=M_n(u[i][1]-k2-n2*Grid[1],n[1]);
            }
        }
        // for Z direction
        for (int  k3 = 0; k3 < Grid[2]; k3++){
            z_direc[i][k3]=M_n(PosIons[i][2]-k3,n[2]);
        }
    }

    fftw_complex *in;   // input variable using standard fftw syntax
    fftw_complex *out;	// output variable
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *Grid[0]*Grid[1]*Grid[2]);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *Grid[0]*Grid[1]*Grid[2]);
    fftw_plan plan;
    plan = fftw_plan_dft_3d(Grid[0],Grid[1],Grid[2], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    #pragma omp for collapse(3) 
    // initializing the "in" vector with zero values
    for (int tx = 0; tx < Grid[0]; tx++){
        for (int ty = 0; ty < Grid[1]; ty++){
            for (int tz = 0; tz < Grid[2]; tz++){
                in[tx * (Grid[2] * Grid[1]) + ty * Grid[2] + tz][REAL] = 0.0;
                // in[tx * (Grid[2] * Grid[1]) + ty * Grid[2] + tz][IMAG] = 0.0;
            }
        }
    }
    #pragma omp for simd
    // Final Q Matrix
    for (int j = 0; j < natoms; j++){
        if (ion_charges[j] == 0)continue;
        for (int tx = 0; tx < Grid[0]; tx++){
            if (x_direc[j][tx] == 0)continue;

            for (int ty = 0; ty < Grid[1]; ty++){
                if (y_direc[j][ty] == 0)continue;

                for (int tz = 0; tz < Grid[2]; tz++){
                    if (z_direc[j][tz] == 0)continue;

                    in[tx * (Grid[2] * Grid[1]) + ty * Grid[2] + tz][0] += ion_charges[j] * x_direc[j][tx] * y_direc[j][ty] * z_direc[j][tz];
                }
            }
        }
    }

    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_cleanup();
    int ii,jj,kk;
    double TwoPi_Grid[] = {2*M_PI/Grid[0],2*M_PI/Grid[1],2*M_PI/Grid[2]};
    double beta_sq = betaa*betaa;
    double pi_sq = M_PI*M_PI;
    double PiSq_BetaSq = pi_sq/beta_sq;
    // #pragma omp parallel for schedule(runtime) reduction(+: reciprocal_energy_i) collapse(3)
    // #pragma omp parallel for schedule(runtime) reduction(+: reciprocal_energy_i) 
    for (int i = -K; i < K+1; i++){
    // for (int i = -Grid[0]/2; i < Grid[0]/2+1; i++){
            double y_part =0;
        for (int j = -K; j< K+1; j++){
        // for (int j = -Grid[1]/2; j< Grid[1]/2+1; j++){
            double z_part =0;
            // for (int k = -K; k < K+1; k++){
            for (int k = -Grid[2]/2; k < Grid[2]/2+1; k++){
                // if(i==0 && j==0 && k==0)continue;
                if(i==0 && j==0)continue;
                if(i<0) ii=Grid[0]+i;
                else ii=i;
                if(j<0) jj=Grid[1]+j;
                else  jj=j;
                if(k<0) kk=Grid[2]+k;
                else  kk=k;
                int temp= ii * (Grid[1] * Grid[2]) + jj * Grid[2] + kk;
                double norm_FQ = out[temp][REAL]*out[temp][REAL]+out[temp][IMAG]*out[temp][IMAG];
                double factor = pow(i/Length[0],2) + pow(j/Length[1],2) + pow((double)k/Grid[2],2);
                // z_part+= norm_FQ * norm(Coeff(TwoPi_Grid[2]*k,n[2])) / (factor*exp(PiSq_BetaSq*factor));
                reciprocal_energy_i+= norm_FQ * norm(Coeff(TwoPi_Grid[0]*i,n[0])*Coeff(TwoPi_Grid[1]*j,n[1])*Coeff(TwoPi_Grid[2]*k,n[2])) / (factor*exp(PiSq_BetaSq*factor));
            }
            // reciprocal_energy_i+=z_part*norm(Coeff(TwoPi_Grid[1]*j,n[1]));
        }
        // reciprocal_energy_i+=y_part*norm(Coeff(TwoPi_Grid[0]*i,n[0]));
    }
    reciprocal_energy_i/=2*M_PI*Length[0]*Length[1]*Grid[2];
    
    // this is the loop for Uo
    for (int  i = 0; i < natoms; i++){
        for (int j = 0; j < natoms; j++){
            reciprocal_energy_o+=ion_charges[i]*ion_charges[j]*F_0(PosIons[i][2]-PosIons[j][2],betaa);
        }
    }
    // return sqrt(M_PI)*reciprocal_energy_o/(Length[0]*Length[1])+reciprocal_energy_i;
    return reciprocal_energy_i;
}