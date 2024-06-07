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

complex<double>Coeff(double v, double w){
    const complex<double> t(0, 1);
    complex<double> bi_mi=exp(t*v*(w-1));
    complex<double> denox;
    for (double f = 0; f < w-1; f++){
        denox+=M_n(f+1,w)*exp(1.0*t*f*v);
    }
    bi_mi/=denox;
    return bi_mi;
}

double reciprocal_fft_integrand(double h, void *params){
    reciprocal_n_params* p = (struct reciprocal_n_params*)params;
    double **PosIons = p->PosIons;
    float *ion_charges = p->ion_charges;
    int natoms = p->natoms;
    double betaa = p->betaa;
    float **box = p->box;
    int K = p->K;
    int Grid = p->Grid;
    // n: order of b-spline interpolation
    int n = p->n;

    // initializing the new variables

    // Structure Factor 
    // complex<double> **StructFact;
    // StructFact = new complex<double> *[Grid];
    // for (int  i = 0; i < Grid; i++){
    //     StructFact[i] = new complex<double> [Grid];
    // }

    // G: Reciprocal Vectors
    // u: the fractional coordinates in x and y directions
    // x_direc, y_direc, z_direc: the cofficients in the x,y and z directions for the Q Matrix
    double G[3][3];
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

    // Length edges of the cell
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    
    int l_max=2;

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
            u[i][j]=Grid*dotProduct(PosIons[i],G[j],3);
        }
    }

    // maximum coordinate along the z axis
    double ZCoord_max = 0;
    // Calculating the cofficients in the x,y and z directions for the Q Matrix
    for (int i = 0; i < natoms; i++){
        // for X direction
        for (int  k1 = 0; k1 < Grid; k1++){
            x_direc[i][k1]=0;
            for (int  n1 = -l_max; n1 < l_max+1; n1++){
                x_direc[i][k1]+=M_n(u[i][0]-k1-n1*K,n);
            }
        }
        // for Y direction
        for (int  k2 = 0; k2 < Grid; k2++){
            y_direc[i][k2]=0;
            for (int  n2 = -l_max; n2 < l_max+1; n2++){
                y_direc[i][k2]+=M_n(u[i][1]-k2-n2*K,n);
            }
        }
        // for Z direction
        for (int  k3 = 0; k3 < Grid; k3++){
            z_direc[i][k3]=M_n(PosIons[i][2]-k3,n);
            ZCoord_max = max(ZCoord_max,PosIons[i][2]); 
        }
    }

    fftw_complex *in;   // input variable using standard fftw syntax
    fftw_complex *out;	// output variable
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *Grid*Grid);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *Grid*Grid);
    fftw_plan plan;
    plan = fftw_plan_dft_2d(Grid,Grid, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    int tz_max = ZCoord_max-n;
    int tz_min = -n;

    for (int x = 0; x < Grid; x++){
        for (int y = 0; y < Grid; y++){
            in[x * Grid + y][0] = 0.0;
            in[x * Grid + y][1] = 0.0;
        }
    }
    // Final Q Matrix after the fourier integral along the z direction
    for (int tz = tz_min; tz <= tz_max; tz++){
        if (natoms == 0) continue;
        for (int j = 0; j < natoms; j++){

            for (int tx = 0; tx < Grid; tx++){
                if (x_direc[j][tx] == 0)continue;

                for (int ty = 0; ty < Grid; ty++){
                    if (y_direc[j][ty] == 0)continue;

                    in[tx * Grid + ty][REAL] += ion_charges[j]*x_direc[j][tx]*y_direc[j][ty]*z_direc[j][tz] * cos(h*tz);
                    in[tx * Grid + ty][IMAG] += ion_charges[j]*x_direc[j][tx]*y_direc[j][ty]*z_direc[j][tz] * sin(h*tz);
                }
            }
        }
    }

    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_cleanup();
 
    long double reciprocal_energy_i=0;
    int ii,jj,kk;
    double deno = 4*betaa*betaa;
    double FourPiPi = 4*M_PI*M_PI;
    double TwoPi_Grid = 2*M_PI/Grid;
    for (int i = -K; i < K+1; i++){
        for (int j = -K; j< K+1; j++){
            if(i==0&&j==0)continue;
            if(i<0) ii=Grid+i;
            else ii=i;
            if(j<0) jj=Grid+j;
            else  jj=j;
            int temp = ii * Grid + jj;
            double factor = FourPiPi * (i*i+j*j) + h*h;
            double norm_FQ = out[temp][REAL]*out[temp][REAL]+out[temp][IMAG]*out[temp][IMAG];
            // cout<<setprecision(18)<<norm_FQ<<"\n";
            // cout<<setprecision(18)<<reciprocal_energy_i<<"\n";
            // cout<<setprecision(30)<<-factor/deno<<"\n";
            // cout<<setprecision(18)<<norm_FQ * exp(-factor/deno) * norm(Coeff(TwoPi_Grid*i,n)*Coeff(TwoPi_Grid*j,n)*Coeff(h,n)) /factor<<"\n";
            // cout<<setprecision(18)<<norm(Coeff(TwoPi_Grid*i,n)*Coeff(TwoPi_Grid*j,n)*Coeff(h,n))<<"\n";
            reciprocal_energy_i+= norm_FQ  * norm(Coeff(TwoPi_Grid*i,n)*Coeff(TwoPi_Grid*j,n)*Coeff(h,n)) / (factor*exp(factor/deno));

        }
    }
    reciprocal_energy_i*=1/(Length[0]*Length[1]);
    return reciprocal_energy_i;
}


double reciprocal_fft(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K, int Grid, int n){
    // this is for Ui
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

    gsl_function F;
    F.function = &reciprocal_fft_integrand; // Set the function to integrate
    reciprocal_n_params params = {PosIons, ion_charges, natoms, betaa, box, K, Grid, n};
    F.params = &params;

    double result, error;
    gsl_integration_qagi(&F, 1e-7, 1e-2, 1000, workspace, &result, &error);
    gsl_integration_workspace_free(workspace); // Free workspace memory
    
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    double reciprocal_energy_o=0;

    // this is the loop for Uo
    for (int  i = 0; i < natoms; i++){
        for (int j = 0; j < natoms; j++){
            reciprocal_energy_o+=ion_charges[i]*ion_charges[j]*F_0(PosIons[i][2]-PosIons[j][2],betaa);
        }
    }
    return result;
    // return sqrt(M_PI)*reciprocal_energy_o/(Length[0]*Length[1])+result;
}