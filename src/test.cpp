#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define REAL 0
#define IMAG 1
const complex<double> t(0.0,1.0);

complex<double> StructureFactor(int mx, int my, double h, double **PosIons, double *ion_charges, int natoms, double **box, int Grid, int n){
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

    complex<double> S;
    double two_pi_mx=2*M_PI*mx,two_pi_my=2*M_PI*my;
    for (int i = 0; i < natoms; i++){
        complex<double> X,Y,Z;
        for (int tx = 0; tx < Grid; tx++){
            X+=x_direc[i][tx]*exp((two_pi_mx*tx)/Grid*t);
        }
        for (int ty = 0; ty < Grid; ty++){
            Y+=y_direc[i][ty]*exp((two_pi_my*ty)/Grid*t);
        }
        for (int  tz = 0; tz < GridZ; tz++){
            Z+=z_direc[i][tz]*exp(h*TZ[tz]*t);
        }
        S+=X*Y*Z*ion_charges[i];
    }
    return S;
}

complex<double> StructureFactor2(int mx, int my, double h, double **PosIons, double *ion_charges, int natoms, double **box, int Grid, int n){
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

    #pragma omp parallel for schedule(runtime)
    // Calculating the fractional coordinates in x and y directions
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 2; j++){
            u[i][j]=Grid*dotProduct(PosIons[i],G[j],3);
        }
    }
    int * TZ = linspace(-nz,(int)Length[2],1);
    int l_max=1;
    #pragma omp parallel for schedule(runtime)
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

    // the fourier integral of z_direc vector for every ith atom
    complex<double>* fz_i_h = new complex<double> [natoms];
    for (int  i = 0; i < natoms; i++){
        for (int  tz = 0; tz < GridZ; tz++){
            fz_i_h[i]+=z_direc[i][tz]*exp(h*TZ[tz]*t);
        }
    }

    // complex<double> * in = new complex<double> [Grid*Grid];
    // for (int i = 0; i < natoms; i++){
    //     for (int tx = 0; tx < Grid; tx++){
    //         if(x_direc[i][tx]==0)continue;
    //         for (int ty = 0; ty < Grid; ty++){   
    //             if(y_direc[i][ty]==0)continue;
    //             in[Grid*tx+ty]+=ion_charges[i]*x_direc[i][tx]*y_direc[i][ty]*fz_i_h[i];
    //         }
    //     }
    // } 

    complex<double> S=0;
    double two_pi_mx=2*M_PI*mx,two_pi_my=2*M_PI*my;
    // for (int tx = 0; tx < Grid; tx++){
    //     for (int ty = 0; ty < Grid; ty++){
    //         S+=in[Grid*tx+ty]*exp((two_pi_my*ty)/Grid*t);
    //     }
    //     S*=exp((two_pi_mx*tx)/Grid*t);
    // }
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

complex<double> StructureFactor3(int mx, int my, double h, double **PosIons, double *ion_charges, int natoms, double **box, int Grid, int n){
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

    #pragma omp parallel for schedule(runtime)
    // Calculating the fractional coordinates in x and y directions
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 2; j++){
            u[i][j]=Grid*dotProduct(PosIons[i],G[j],3);
        }
    }
    int * TZ = linspace(-nz,(int)Length[2],1);
    int l_max=1;
    #pragma omp parallel for schedule(runtime)
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
    fftw_destroy_plan(plan);
    fftw_cleanup();
    if(mx<0)mx=Grid+mx;
    if(my<0)my=Grid+my;
    complex<double> S=0;
    S = out[Grid*mx+my][0]+ t*out[Grid*mx+my][1];
    return S;
}