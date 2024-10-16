// this code implements the reciprocal (k!=0) energy using the bspline method with 2D Fourier Transform and 1D Fourier Integral
#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define REAL 0
#define IMAG 1
const complex<double> t(0.0,1.0);

#define ENABLE_OMP 1

struct reciprocal_n_params {
    double** PosIons;
    double* ion_charges;
    int natoms;
    double betaa;
    double** box;
    int K;
    int Grid;
    int n;
    double* Length;
    double **G;
    double **x_direc, **y_direc, **z_direc;
    int * TZ ;
    int GridZ;
};

complex<double> func2(int mx, int my, int Grid, double **x_direc, double **y_direc, double *ion_charges, int natoms, complex<double>* fz_i_h){
    complex<double> S;
    double two_pi_mx=2*M_PI*mx,two_pi_my=2*M_PI*my;
    for(int tx = 0; tx < Grid; tx++){
        complex<double> Sx=0;
        for (int ty = 0; ty < Grid; ty++){
            complex<double> Sy=0;
            for (int i = 0; i < natoms; i++){
                Sy+=ion_charges[i]*x_direc[i][tx]*y_direc[i][ty]*fz_i_h[i];
            }
            Sx+=Sy*exp((two_pi_my*ty)/Grid*t);
        }
        S+=Sx*exp((two_pi_mx*tx)/Grid*t);
    }

    return S;
}

double reciprocal_ft_integrand(double h, void *params){
    reciprocal_n_params* p = (struct reciprocal_n_params*)params;
    double **PosIons = p->PosIons;
    double *ion_charges = p->ion_charges;
    int natoms = p->natoms;
    double betaa = p->betaa;
    double **box = p->box;
    int K = p->K;
    int Grid = p->Grid;
    // n: order of b-spline interpolation
    int n = p->n;
    double *Length = p->Length;
    auto G = p->G;
    auto x_direc=p->x_direc,y_direc=p->y_direc,z_direc=p->z_direc;
    int *TZ=p->TZ;
    int GridZ = p->GridZ;

    long double reciprocal_energy_i=0;
    double deno = 4*betaa*betaa;
    double FourPiPi = 4*M_PI*M_PI;
    double TwoPi_Grid = 2*M_PI/Grid;

    // the fourier integral of z_direc vector for every ith atom
    complex<double>* fz_i_h = new complex<double> [natoms];

    for (int  i = 0; i < natoms; i++){
        for (int  tz = 0; tz < GridZ; tz++){
            fz_i_h[i]+=z_direc[i][tz]*exp(h*TZ[tz]*t);
        }
    }

    fftw_complex *in;
    fftw_complex *out;

    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Grid*Grid);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Grid*Grid);

    fftw_plan plan;
    plan = fftw_plan_dft_2d(Grid, Grid, in ,out, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int i = 0; i < natoms; i++){
        for (int tx = 0; tx < Grid; tx++){
            if(x_direc[i][tx]==0)continue;
            for (int ty = 0; ty < Grid; ty++){   
                if(y_direc[i][ty]==0)continue;
                in[Grid*tx+ty][0]+=ion_charges[i]*x_direc[i][tx]*y_direc[i][ty]*fz_i_h[i].real();
                in[Grid*tx+ty][1]+=ion_charges[i]*x_direc[i][tx]*y_direc[i][ty]*fz_i_h[i].imag();
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
            if(i<0)ii=Grid+i;
            else ii=i;
            if(j<0)jj=Grid+j;
            else jj=j;
            int temp = Grid*ii+jj;
            double factor = FourPiPi * (i*i/(box[0][0]*box[0][0])+j*j/(box[1][1]*box[1][1])) + h*h;

            // double norm_F = norm(func2(i, j, Grid, x_direc,y_direc,ion_charges,natoms,fz_i_h));
            double norm_FQ = norm(out[temp][0] + t*out[temp][1]);
            // if(abs(norm_F-norm_FQ)>0.01){
            //     cout<<norm_FQ<<" "<<norm_F<<" "<<out[temp][0] + t*out[temp][1]<<" "<<i<<" "<<j<<" "<<h<<"\n";
            // }
            reciprocal_energy_i+= norm_FQ  * norm(Coeff(TwoPi_Grid*i,n)*Coeff(TwoPi_Grid*j,n)*Coeff(h,8)) / (factor*exp(factor/deno));
        }
    }
    return reciprocal_energy_i;
}

double reciprocal_fft(double **PosIons, double *ion_charges, int natoms, double betaa, double **box, int K, int Grid, int n){
    // this is for Ui
    // Edge lengths of the cell
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    // Volume Calculations
    double A[3];
    double C[3]={box[2][0],box[2][1],box[2][2]};
    crossProduct(box[0],box[1],A);
    double volume = dotProduct(A,C,3);
    int nz = 8;
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
    int GridZ = Length[2]+nz+1;
    double **u,**x_direc, **y_direc, **z_direc;
    u= new double * [natoms];
    x_direc = new double * [natoms];
    y_direc = new double * [natoms];
    z_direc = new double * [natoms]; 

    for (int  i = 0; i < natoms; i++){
        u[i] = new double  [2]; // We only need these in x and y direction 
        x_direc[i] = new double  [Grid];
        y_direc[i] = new double  [Grid];
        z_direc[i] = new double  [GridZ]; // tz varies from -n to Zmax(Lz)
    }

    // Calculating the fractional coordinates in x and y directions
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 2; j++){
            u[i][j]=Grid*dotProduct(PosIons[i],G[j],3);
        }
    }
    int * TZ = linspace(-nz,(int)Length[2],1);
    int l_max=1;

    // Calculating the cofficients in the x,y and z directions for the Q Matrix
    for (int i = 0; i < natoms; i++){
        // for X direction
        for (int  tx = 0; tx < Grid; tx++){
            x_direc[i][tx]=0;
            for (int  lx = -l_max; lx < l_max+1; lx++){
                x_direc[i][tx]+=M_n(u[i][0]-tx-lx*Grid,n);
            }
        }
        // for Y direction
        for (int  ty = 0; ty < Grid; ty++){
            y_direc[i][ty]=0;
            for (int  ly = -l_max; ly < l_max+1; ly++){
                y_direc[i][ty]+=M_n(u[i][1]-ty-ly*Grid,n);
            }
        }
        // for Z direction
        for (int  tz = 0; tz < GridZ; tz++){
            z_direc[i][tz]=M_n(PosIons[i][2]-TZ[tz],nz);
        }
    }
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(200);
    // gsl_integration_romberg_workspace *workspace = gsl_integration_romberg_alloc(12);

    gsl_function F;
    F.function = &reciprocal_ft_integrand; // Set the function to integrate
    reciprocal_n_params params = {PosIons, ion_charges, natoms, betaa, box, K, Grid, n, Length, G,x_direc,y_direc,z_direc, TZ, GridZ};
    F.params = &params;
    double result, error;
    // size_t size = 12;
    // gsl_integration_romberg(&F,-10,10, 1e-7, 1e-2, &result, &size, workspace);
    // gsl_integration_qag(&F, -5, 5, 1e-4, 1e-2, 200, GSL_INTEG_GAUSS15, workspace, &result, &error);
    gsl_integration_qagi(&F, 1e-4, 1e-2, 200, workspace, &result, &error);
    // gsl_integration_romberg_free(workspace); // Free workspace memory
    gsl_integration_workspace_free(workspace); // Free workspace memory
    
    result*=1/(Length[0]*Length[1]);
    return result;
}