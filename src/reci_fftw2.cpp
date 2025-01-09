// this code implements the reciprocal (k!=0) energy using the bspline method with 2D Fourier Transform and 1D Fourier Integral
#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define REAL 0
#define IMAG 1
// #define ENABLE_OMP 1
const complex<double> t(0.0,1.0);
double FourPiPi = 4*M_PI*M_PI;

struct reciprocal_n_params {
    double* ion_charges;
    int natoms;
    double betaa;
    int K;
    int *Grid;
    int *n;
    double **G;
    double **x_direc, **y_direc, **z_direc;
    double * TZ ;
    int GridZ;
};

double reciprocal_ft_integrand(double h, void *params){
    reciprocal_n_params* p = (struct reciprocal_n_params*)params;
    double *ion_charges = p->ion_charges;
    int natoms = p->natoms;
    double betaa = p->betaa;
    int K = p->K;
    int *Grid = p->Grid;
    // n: order of b-spline interpolation
    int *n = p->n;
    auto G = p->G;
    auto x_direc=p->x_direc,y_direc=p->y_direc,z_direc=p->z_direc;
    double *TZ=p->TZ;
    int GridZ = p->GridZ;

    long double reciprocal_energy_i=0;
    double deno = 4*betaa*betaa;
    double TwoPi_Gridx = 2*M_PI/Grid[0];
    double TwoPi_Gridy = 2*M_PI/Grid[1];

    // the fourier integral of z_direc vector for every ith atom
    complex<double>* fz_i_h = new complex<double> [natoms];
    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime)
    #endif
    for (int  i = 0; i < natoms; i++){
        for (int  tz = 0; tz < GridZ; tz++){
            fz_i_h[i]+=z_direc[i][tz]*exp(h*TZ[tz]*t);
        }
    }

    fftw_complex *in;
    fftw_complex *out;

    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Grid[0]*Grid[1]);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Grid[0]*Grid[1]);

    fftw_plan plan;
    plan = fftw_plan_dft_2d(Grid[0], Grid[1], in ,out, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int i = 0; i < natoms; i++){
        for (int tx = 0; tx < Grid[0]; tx++){
            if(x_direc[i][tx]==0)continue;
            for (int ty = 0; ty < Grid[1]; ty++){   
                if(y_direc[i][ty]==0)continue;
                in[Grid[1]*tx+ty][0]+=ion_charges[i]*x_direc[i][tx]*y_direc[i][ty]*fz_i_h[i].real();
                in[Grid[1]*tx+ty][1]+=ion_charges[i]*x_direc[i][tx]*y_direc[i][ty]*fz_i_h[i].imag();
            }
        }
    }

    fftw_execute(plan);
    // fftw_destroy_plan(plan);
    // fftw_cleanup();

    int ii,jj;
    for (int i = -K; i < K+1; i++){
        for (int j = -K; j< K+1; j++){
            if(i==0&&j==0)continue;
            if(i<0)ii=Grid[0]+i;
            else ii=i;
            if(j<0)jj=Grid[1]+j;
            else jj=j;
            int temp = Grid[1]*ii+jj;
            double factor = FourPiPi * (i*i*G[0][0]*G[0][0]+j*j*G[1][1]*G[1][1]) + h*h;
            double norm_FQ = norm(out[temp][0] + t*out[temp][1]);
            reciprocal_energy_i+= norm_FQ  * norm(Coeff(TwoPi_Gridx*i,n[0])*Coeff(TwoPi_Gridy*j,n[1])*Coeff(h,8)) / (factor*exp(factor/deno));
        }
    }
    return reciprocal_energy_i;
}

// Main Function to Calculate the Reciprocal Energy (k!=0)
double reciprocal_fft(double **PosIons, double *ion_charges, int natoms, double betaa, double **box, int K, int *Grid, int *n){

    #if defined ENABLE_OMP
        omp_set_num_threads(thread::hardware_concurrency());
    #endif
    
    // Edge lengths of the cell
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};

    // Volume Calculations
    double A[3];
    double C[3]={box[2][0],box[2][1],box[2][2]};
    crossProduct(box[0],box[1],A);
    double volume = dotProduct(A,C,3);

    // Calculating the reciprocal vectors
    double **G;
    G= new double * [3];
    for (int i = 0; i < 3; i++){
        G[i] = new double  [3];
    }
    crossProduct(box[1],box[2],G[0]);
    crossProduct(box[2],box[0],G[1]);
    crossProduct(box[0],box[1],G[2]);
    for (int x = 0; x < 3; x++)
        for (int q = 0; q < 3; q++)
            G[x][q] /= volume;

    // initializing the new variables
    // u: the fractional coordinates in x and y directions
    // x_direc, y_direc, z_direc: the cofficients in the x,y and z directions for the Q Matrix
    int GridZ = Length[2]+n[2]+1;
    double **u,**x_direc, **y_direc, **z_direc;
    u= new double * [natoms];
    x_direc = new double * [natoms];
    y_direc = new double * [natoms];
    z_direc = new double * [natoms]; 
    for (int  i = 0; i < natoms; i++){
        u[i] = new double  [2]; // We only need these in x and y direction 
        x_direc[i] = new double  [Grid[0]];
        y_direc[i] = new double  [Grid[1]];
        z_direc[i] = new double  [GridZ]; // tz varies from -n to Zmax(Lz)
    }

    // Calculating the fractional coordinates in x and y directions
    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime)
    #endif
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 2; j++){
            u[i][j]=Grid[j]*dotProduct(PosIons[i],G[j],3);
        }
    }
    double * TZ = linspace(-n[2],(int)Length[2],1);
    int l_max=1;

    // Calculating the cofficients in the x,y and z directions for the Q Matrix
    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime)
    #endif
    for (int i = 0; i < natoms; i++){
        // for X direction
        for (int  tx = 0; tx < Grid[0]; tx++){
            x_direc[i][tx]=0;
            for (int  lx = -l_max; lx < l_max+1; lx++){
                x_direc[i][tx]+=M_n(u[i][0]-tx-lx*Grid[0],n[0]);
            }
        }
        // for Y direction
        for (int  ty = 0; ty < Grid[1]; ty++){
            y_direc[i][ty]=0;
            for (int  ly = -l_max; ly < l_max+1; ly++){
                y_direc[i][ty]+=M_n(u[i][1]-ty-ly*Grid[1],n[1]);
            }
        }
        // for Z direction
        for (int  tz = 0; tz < GridZ; tz++){
            z_direc[i][tz]=M_n(PosIons[i][2]-TZ[tz],n[2]);
        }
    }
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(200);
    gsl_function F;
    F.function = &reciprocal_ft_integrand; // Set the function to integrate
    reciprocal_n_params params = {ion_charges, natoms, betaa, K, Grid, n, G, x_direc, y_direc, z_direc, TZ, GridZ};
    F.params = &params;
    double result, error;
    gsl_integration_qagi(&F, 1e-4, 1e-2, 200, workspace, &result, &error);
    gsl_integration_workspace_free(workspace); // Free workspace memory
    // gsl_integration_romberg_workspace *workspace = gsl_integration_romberg_alloc(12);
    // size_t size = 12;
    // gsl_integration_romberg(&F,-10,10, 1e-7, 1e-2, &result, &size, workspace);
    // gsl_integration_romberg_free(workspace); // Free workspace memory
    
    // gsl_integration_qag(&F, -5, 5, 1e-4, 1e-2, 200, GSL_INTEG_GAUSS15, workspace, &result, &error);

    result*=1/(Length[0]*Length[1]);
    return result;
}