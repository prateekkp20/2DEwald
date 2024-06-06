#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define REAL 0
#define IMAG 1

struct reciprocal_n_params {
  double** PosIons;
  float* ion_charges;
  int natoms;
  double betaa;
  float** box;
  int K;
  int Grid;
  int n;
};

double M_n(double u, int n){
    if(n<2)return 0;
    else if(n==2){
        if(u<0 || u>n) return 0;
        else{
            return 1-abs(u-1);
        }
    }
    else{
        if(u<0 || u>n) return 0;
        else{
            return (u*M_n(u,n-1)/(n-1))+((n-u)*M_n(u-1,n-1)/(n-1));
        }
    }
}

// double reciprocal_fftw_integrand(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K, int Grid, int n){
double reciprocal_fftw_integrand(double h, void *params){
    reciprocal_n_params* p = (struct reciprocal_n_params*)params;
    double **PosIons = p->PosIons;
    float *ion_charges = p->ion_charges;
    int natoms = p->natoms;
    double betaa = p->betaa;
    float **box = p->box;
    int K = p->K;
    int Grid = p->Grid;
    int n = p->n;

    // n: order of b-spline interpolation
    // initializing the new variables
    // Structure Factor 
    complex<double> **StructFact;
    StructFact = new complex<double> *[Grid];
    for (int  i = 0; i < Grid; i++){
        StructFact[i] = new complex<double> [Grid];
    }

    for (int  i = 0; i < Grid; i++){
        StructFact[i] = new complex<double> *[Grid];
        for (int j = 0; j < Grid; j++){
            StructFact[i][j] = new complex<double> [Grid];
        }
    }
    cout<<StructFact[1][2][2]<<"\n";

    double G[3][3], m[3];
    double **u,**x_direc, **y_direc, **z_direc;
    u= new double * [natoms];
    x_direc = new double * [natoms];
    y_direc = new double * [natoms];
    z_direc = new double * [natoms];

    for (int  i = 0; i < natoms; i++){
        u[i] = new double  [2]; // We only need these in x and y direction 
        x_direc[i] = new double  [Grid];
        y_direc[i] = new double  [Grid];
        z_direc[i] = new double  [Grid];
    }

    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    int n_max=2;

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

    // Calculating the fractional coordinates of x and y
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 2; j++){
            u[i][j]=Grid*dotProduct(PosIons[i],G[j],3);
        }
    }
    
    // Calculating the cofficients in the x,y and z directions for the Q Matrix
    for (int i = 0; i < natoms; i++){
        // for X direction
        for (int  k1 = 0; k1 < Grid; k1++){
            x_direc[i][k1]=0;
            for (int  n1 = -n_max; n1 < n_max+1; n1++){
                x_direc[i][k1]+=M_n(u[i][0]-k1-n1*K,n);
            }
        }
        // for Y direction
        for (int  k2 = 0; k2 < Grid; k2++){
            y_direc[i][k2]=0;
            for (int  n2 = -n_max; n2 < n_max+1; n2++){
                y_direc[i][k2]+=M_n(u[i][1]-k2-n2*K,n);
            }
        }
        // for Z direction
        for (int  k3 = 0; k3 < Grid; k3++){
            z_direc[i][k3]=M_n(PosIons[i][2]-k3,n);
        }
    }

    fftw_complex *in;   // input variable using standard fftw syntax
    fftw_complex *out;	// output variable
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *Grid*Grid);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *Grid*Grid);
    fftw_plan p;
    p = fftw_plan_dft_2d(K,K, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int i = 0; i < Grid; i++){
        
    }
    


    double reciprocal_energy_o=0;
    // this is the loop for Uo
    #pragma omp parallel for simd schedule(runtime) reduction(+: reciprocal_energy_o) collapse(2)
    for (int  i = 0; i < natoms; i++){
        for (int j = 0; j < natoms; j++){
            reciprocal_energy_o+=ion_charges[i]*ion_charges[j]*F_0(PosIons[i][2]-PosIons[j][2],betaa);
        }
    }
    return sqrt(M_PI)*reciprocal_energy_o/(Length[0]*Length[1]);
}