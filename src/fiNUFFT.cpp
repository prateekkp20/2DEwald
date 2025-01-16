#include "libinclude.h"
#include "const.h"
#include "fundec.h"

#define REAL 0
#define IMAG 1
const complex<double> t(0.0,1.0);

double without(double **PosIons, int natoms, double **box, int *n){
    omp_set_num_threads(thread::hardware_concurrency());
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    int GridZ = Length[2]+n[2]+1;
    double **z_direc;
    z_direc = new double * [natoms]; 
    for (int  i = 0; i < natoms; i++){
        z_direc[i] = new double  [GridZ]; // tz varies from -n to Zmax(Lz)
    }

    double * TZ = linspace(-n[2],(int)Length[2],1);
    // Calculating the cofficients in the x,y and z directions for the Q Matrix
    #pragma omp parallel for
    for (int i = 0; i < natoms; i++){
        // for Z direction
        for (int  tz = 0; tz < GridZ; tz++){
            z_direc[i][tz]=M_n(PosIons[i][2]-TZ[tz],n[2]);
        }
    }

    vector<double> H;
    double start = -5;
    for (int i = 0; i < GridZ; i++){
        H.push_back(start+0.5);
    }

    // the fourier integral of z_direc vector for every ith atom
    complex<double> **fz_i_h;
    fz_i_h = new complex<double> * [natoms];
    #pragma omp parallel for
    for (int  i = 0; i < natoms; i++){
        fz_i_h[i] = new complex<double> [H.size()];
        int count = 0;
        for(double itr: H){
            for (int  tz = 0; tz < GridZ; tz++){
                fz_i_h[i][count]+=z_direc[i][tz]*exp(itr*TZ[tz]*t);
            }
            count++;
        }
    }
    return 0;
}

double with(double **PosIons, int natoms, double **box, int *n){
    double Length[3]={sqrt(dotProduct(box[0],box[0],3)),sqrt(dotProduct(box[1],box[1],3)),sqrt(dotProduct(box[2],box[2],3))};
    int GridZ = Length[2]+n[2]+1;
    double **z_direc;
    z_direc = new double * [natoms]; 
    for (int  i = 0; i < natoms; i++){
        z_direc[i] = new double  [GridZ]; // tz varies from -n to Zmax(Lz)
    }

    double * TZ = linspace(-n[2],(int)Length[2],1);
    // Calculating the cofficients in the x,y and z directions for the Q Matrix
    for (int i = 0; i < natoms; i++){
        // for Z direction
        for (int  tz = 0; tz < GridZ; tz++){
            z_direc[i][tz]=M_n(PosIons[i][2]-TZ[tz],n[2]);
        }
    }
    int type = 3, dim = 1;
    
    return 0;
}