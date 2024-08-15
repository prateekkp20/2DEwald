//***********************************************
// function to calculate distances according to minimum image convention for a 2D symmetric system
//***********************************************

#include "libinclude.h" 

double dist(double **PosIons, int atom1, int atom2, float **box){

	double Dx, Dy, Dz;
	double Dx1, Dy1, Dz1;
    double distatoms;
	
	Dx = PosIons[atom1][0] - PosIons[atom2][0];
	Dy = PosIons[atom1][1] - PosIons[atom2][1];
	Dz = PosIons[atom1][2] - PosIons[atom2][2];
	
    Dx1 = Dx - box[0][0]*ceil(Dx/box[0][0]-0.5) - box[1][0]*ceil(Dy/box[1][1]-0.5);
	Dy1 = Dy - box[0][1]*ceil(Dx/box[0][0]-0.5) - box[1][1]*ceil(Dy/box[1][1]-0.5);
	Dz1 = Dz - box[0][2]*ceil(Dx/box[0][0]-0.5) - box[1][2]*ceil(Dy/box[1][1]-0.5);
	
	distatoms = sqrt(pow(Dx1,2) + pow(Dy1,2) + pow(Dz1,2));
	return distatoms;
	
}

// this function returns the modR, deltaZ and the deltaZ^2 in a three element vector 
// void distWithZ(double **PosIons, int atom1, int atom2, float **box, vector<double> &out){
void distWithZ(double **PosIons, int atom1, int atom2, float **box, double &modR, double &Z, double &Z2){
	double Dx, Dy, Dz;
	double Dx1, Dy1;
	
	Dx = PosIons[atom1][0] - PosIons[atom2][0];
	Dy = PosIons[atom1][1] - PosIons[atom2][1];
	Dz = PosIons[atom1][2] - PosIons[atom2][2];
	
    Dx1 = Dx - box[0][0]*ceil(Dx/box[0][0]-0.5) - box[1][0]*ceil(Dy/box[1][1]-0.5);
	Dy1 = Dy - box[0][1]*ceil(Dx/box[0][0]-0.5) - box[1][1]*ceil(Dy/box[1][1]-0.5);
	Z = Dz - box[0][2]*ceil(Dx/box[0][0]-0.5) - box[1][2]*ceil(Dy/box[1][1]-0.5);

	Z2 = pow(Z,2);
	modR = sqrt(pow(Dx1,2) + pow(Dy1,2) + Z2);

	// vector<double> out(3,0);
	// modR=distatoms;
	// return {distatoms,Dz1,Dz2};
	
}
