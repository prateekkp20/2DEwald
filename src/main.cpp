#include "libinclude.h" 
#include "potdec.h"
#include "mathdec.h"
#include "fundec.h"
#include "const.h"
#include "fftw3.h"
#include "complex"
#include "header.h"

double *ExpFactor;
double G[3][3];
double volume;
complex<double> *CoeffX,*CoeffY,*CoeffZ;

int main(int argc, char **argv){

	string filename;
	string posfile;
	string garbage, garbage1;
	char garbage2[200];
	int tmp;
	int i, j, MC, MD;
	int  nMD, md_print;
	float Temp;
	int  MDtime, MDeq;
	float dt;
	int genpairlist;
	int randnum;
	int hyb, NMCMD, rdvel;
	float tmp2;
	double unitzer=(AVG*COL*CHARELEC*CHARELEC/rtoa)/1000/4.184;

	//--------------------get the name of the file for input variables, default is input.in-----------------------//
	if(argc <2)
		filename = "input.in";
	else
		filename = argv[1];

	ifstream InputIn(filename.c_str(),ios::in);
	if(!InputIn){
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}

	InputIn>>garbage>>garbage;
	getline(InputIn, garbage);

	// cout<<"Project Name: "<<garbage<<endl;

	InputIn>>garbage>>garbage;
	InputIn>>posfile;

	InputIn>>garbage>>garbage;
	InputIn>>genpairlist;

	InputIn>>garbage>>garbage;
	InputIn>>Temp;

	// cout<<"Temperature of the system:\t\t"<<Temp<<" K"<<endl;

	InputIn>>garbage>>garbage;
	InputIn>>hyb;

	InputIn>>garbage>>garbage;
	InputIn>>NMCMD;

	InputIn>>garbage>>garbage;
	InputIn>>rdvel;

	bool comdens;

	InputIn>>garbage>>garbage;
	InputIn>>comdens;

	// cout<<"Compute Density: "<<comdens<<"\n";

	double lim1, lim2, dx;

	InputIn>>garbage>>garbage;
	InputIn>>lim1>>lim2>>dx;

	// cout<<"Range of x-values considered for computing densities: "<<lim1<<"\t"<<lim2<<endl;

	//Parameters for simulated annealing		
	bool simann;

	InputIn>>garbage>>garbage;
	InputIn>>simann;

	float TempB, TempE;

	InputIn>>garbage>>garbage;
	InputIn>>TempB>>TempE;

	int rampsize, rampstep;

	InputIn>>garbage>>garbage;
	InputIn>>rampsize>>rampstep;

	int ncycle;

	InputIn>>garbage>>garbage;
	InputIn>>ncycle;

	if(simann==1){	
		cout<<"Perform Simulated Annealing"<<endl;
		cout<<"Start Temperature: "<<TempB<<" K"<<endl;
		cout<<"End Temperature: "<<TempE<<" K"<<endl;
		cout<<"Delta T:"<<rampsize<<" K"<<endl;
		cout<<"Increase temperature after: "<<rampstep<<" MD steps"<<endl;
		cout<<"ncycles: "<<ncycle<<endl;
	}

	bool minimize;

	InputIn>>garbage>>garbage;
	InputIn>>minimize;                      //whether to minimize or not
	if(minimize)cout<<"****Perform minimization****"<<endl;

	InputIn.close();

	//--------------------------------------- ewald input file ---------------------------------------//

	ifstream EWALDIn("ewald.in",ios::in);
	if(!EWALDIn){
		cerr << "File EWALDIn could not be opened" <<endl;
		exit(1);
	}

	// convergence limits for the reci sum
	// grid size for x and y direction
	// order of bspline interpolation in the x, y and z direction
	int Kvec[3], grid[3], order[3];
	for (int i = 0; i < 3; i++){
		EWALDIn>>garbage>>garbage;
		EWALDIn>>Kvec[i];
	}

	for (int i = 0; i < 3; i++){
		EWALDIn>>garbage>>garbage;
		EWALDIn>>grid[i];
	}

	for (int i = 0; i < 3; i++){
		EWALDIn>>garbage>>garbage;
		EWALDIn>>order[i];
	}                    

	EWALDIn.close();

	// cout<<"Initial positions read from the file: "<<posfile<<endl;

	//--------------------------------------- get positions from the POSCAR file ---------------------------------------//

	ifstream PosIn(posfile.c_str(),ios::in);
	if(!PosIn){
		cerr << "Posfile could not be opened" <<endl;
		exit(1);
	}

	getline(PosIn, garbage);
	istringstream StrStream(garbage);

	// Counting through the types of the atoms present in the cell and storing in the n_atomtype variable
	int n_atomtype=0;
	while(StrStream){
		getline(StrStream, garbage1, ' ');
		if(garbage1.compare("") != 0)
			n_atomtype = n_atomtype + 1;
	}
	
	string *atomtype;
	atomtype=new string [n_atomtype];

	istringstream strstream(garbage);

	tmp = 0;

	while(strstream){
		getline(strstream, garbage1, ' ');
		if(garbage1.compare("") != 0){
			atomtype[tmp]=garbage1;
			tmp = tmp + 1;
		}
	}
	
	int *natoms_type, natoms;  //natoms_type: number of atoms for each type, natoms: total atoms in the unit cell

	natoms_type=new int [n_atomtype];
	natoms=0;

	//Creating 3x3 boxcell array
	double **boxcell;
	boxcell=new double * [3];
	for(i = 0; i<3;i++){
		boxcell[i]=new double [3];
	}

	getline(PosIn, garbage); //read overall scaling factor

	//get box cell vectors
	for(i = 0; i<3;i++){
		for(j=0;j<3;j++){
			PosIn>>boxcell[i][j];
		}
	}
	getline(PosIn, garbage);//read the atoms again
	getline(PosIn, garbage1);
	
	//get number of atoms for each type
	for(i=0;i<n_atomtype;i++){
		PosIn>>natoms_type[i];
		natoms=natoms+natoms_type[i];
	}

	getline(PosIn, garbage);
	getline(PosIn, garbage);

	double **PosIons, *PosIons2, *ForceIons, *vel;
	int *fixatoms;
	double *mass;

	mass=new double [natoms];
	PosIons=new double * [natoms];
	PosIons2=new double [natoms*3];
	ForceIons=new double [natoms*3];
	vel=new double [natoms*3];
	fixatoms=new int [natoms*3];

	for(i=0;i<natoms;i++){
		PosIons[i]=new double [3];
	}

	//getting positions of each atom
	for(i=0;i<natoms;i++){
		PosIn>>PosIons[i][0]>>PosIons[i][1]>>PosIons[i][2];
		PosIons2[3*i]=PosIons[i][0];
		PosIons2[3*i+1]=PosIons[i][1];
		PosIons2[3*i+2]=PosIons[i][2];
	}

	PosIn.close();

	//--------------------------------------- charge input file ---------------------------------------//

	ifstream CHARGEIn("charge.in",ios::in);
	if(!CHARGEIn){
		cerr << "File CHARGEIn could not be opened" <<endl;
		exit(1);
	}

	double *chg;
	chg=new double [n_atomtype];

	for(int i=0; i<n_atomtype; i++){
		CHARGEIn>>chg[i];
	}

	CHARGEIn.close();

	/*Checking Charge Neutrality for the system*/
	double total_charge = 0;
	for (int i = 0; i < n_atomtype; i++){
		total_charge += chg[i]* natoms_type[i];
	}
	if(total_charge){
		cout<<"Error: System is not charge neutral"<<endl;
		return 0;
	}
	
	//-------------------------creating the charge array for each atom present in the unit cell ------------------------------//
	double *ion_charges;
	double *charge_prod;
	ion_charges=new double [natoms];
	charge_prod=new (nothrow) double [natoms*(natoms-1)/2];
	if (!charge_prod){
    	cout << "Memory allocation failed\n";
	}

	int c=0;
	for (int i = 0; i < n_atomtype; i++){
		for (int j = 0; j < natoms_type[i]; j++){
			ion_charges[c]=chg[i];
			c++;
		}		
	}
	for (int i = 1; i < natoms; i++){
		for (int j = 0; j < i; j++){
			charge_prod[i*(i-1)/2+j] = ion_charges[i]*ion_charges[j];
		}
	}

	// print_carcoor(PosIons2, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "CONTCAR");
	// print_coor(PosIons2, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "COOR");
	// print_lammps_input_file(PosIons2, chg, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "data.lmp");

	/*Check for orthogonality of sides*/
	/*Ewald method only works for unit cell with orthogonal sides*/
	if(dotProduct(boxcell[1],boxcell[0],3) || dotProduct(boxcell[2],boxcell[0],3) || dotProduct(boxcell[1],boxcell[2],3)){
		cout<<"Error: Unit Cell with Non Orthogonal Sides"<<endl;
		return 0;
	}

	double Lmin=min(boxcell[0][0],min(boxcell[1][1],boxcell[2][2]));
	double beta=5.42/Lmin;
	double cutoff = Lmin/2;

	// Volume Calculations
    double A[3];
    double C[3]={boxcell[2][0],boxcell[2][1],boxcell[2][2]};
    crossProduct(boxcell[0],boxcell[1],A);
    double volume = dotProduct(A,C,3);
	// Calculating the reciprocal vectors
    crossProduct(boxcell[1],boxcell[2],G[0]);
    crossProduct(boxcell[2],boxcell[0],G[1]);
    crossProduct(boxcell[0],boxcell[1],G[2]);
    for (int x = 0; x < 3; x++)
        for (int q = 0; q < 3; q++)
            G[x][q] /= volume;

	/*Useful configuration independent Computations for the reciprocal space summation*/
	/*Cofficients Bi[mi] of the bspline interpolation*/ //Refer to Essmann et al.
	// The negative indices are stored from the end of the array
	CoeffX = new complex<double> [2*Kvec[0]+1];
    CoeffY = new complex<double> [2*Kvec[1]+1];
    CoeffZ = new complex<double> [2*Kvec[2]+1];
	double TwoPi_Gridx = 2*M_PI/grid[0];
    double TwoPi_Gridy = 2*M_PI/grid[1];
    double TwoPi_Gridz = 2*M_PI/grid[2];
	for (int i = -Kvec[0]; i < Kvec[0]+1; i++){
        int ic;
        if(i<0) {ic=(2*Kvec[0]+1)+i;}
        else {ic=i;}
        CoeffX[ic] = Coeff(TwoPi_Gridx*i,order[0]);
    }
	for (int i = -Kvec[1]; i < Kvec[1]+1; i++){
        int ic;
        if(i<0) {ic=(2*Kvec[1]+1)+i;}
        else {ic=i;}
        CoeffY[ic] = Coeff(TwoPi_Gridy*i,order[1]);
    }
	for (int i = -Kvec[2]; i < Kvec[2]+1; i++){
        int ic;
        if(i<0) {ic=(2*Kvec[2]+1)+i;}
        else {ic=i;}
        CoeffZ[ic] = Coeff(TwoPi_Gridz*i,order[2]);
    }

	/* B(m1,m2,m3)*Exp(-|G|)/|G| term in the reciprocal loop*/
    double L1 = boxcell[0][0];
    double L2 = boxcell[1][1];
    double L3 = boxcell[2][2];
	ExpFactor = new double [(2*Kvec[0]+1)*(2*Kvec[1]+1)*(2*Kvec[2]+1)];
	#pragma omp parallel for schedule(runtime) collapse(3)
	for (int i = -Kvec[0]; i < Kvec[0]+1; i++){
        for (int j = -Kvec[1]; j< Kvec[1]+1; j++){
            for (int k = -Kvec[2]; k < Kvec[2]+1; k++){
				if(i==0&&j==0&&k==0)continue;
				int ii,jj,kk;
				if(i<0) ii=(2*Kvec[0]+1)+i;
                else ii=i;
                if(j<0) jj=(2*Kvec[1]+1)+j;
                else  jj=j;
                if(k<0) kk=(2*Kvec[2]+1)+k;
                else  kk=k;
				int temp=ii * ((2*Kvec[2]+1) * (2*Kvec[1]+1)) + jj * (2*Kvec[2]+1) + kk;
				ExpFactor[temp] = norm(CoeffX[ii]*CoeffY[jj]*CoeffZ[kk])*constantterm(i,j,k,L1,L2,L3,beta,600);
			}
		}
	}

	/*Self Energy*/
	double selfenergy=self(n_atomtype, natoms_type, chg, beta)*unitzer;
	cout<<fixed<<setprecision(5)<<"Self Energy: "<<selfenergy<<" Kcal/mol"<<"\n\n";

	/*Reciprocal Energy (k!=0)*/
	chrono::time_point<std::chrono::system_clock> start1, end1;
	start1 = chrono::system_clock::now();
	double recienergy=reciprocal_n2(PosIons, charge_prod, ion_charges, natoms, beta, boxcell, Kvec)*unitzer;
	cout<<fixed<<setprecision(15)<<"Reciprocal Energy: "<<recienergy<<" Kcal/mol"<<"\n";
	end1 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds1 = end1- start1;
    time_t end_time1 = std::chrono::system_clock::to_time_t(end1);
	cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds1.count() << " sec\n\n";

	/*Real Energy*/
	chrono::time_point<std::chrono::system_clock> start2, end2;
	start2 = chrono::system_clock::now();
	double realenergy=real(PosIons2, charge_prod, natoms, beta, boxcell,cutoff)*unitzer;
	cout<<fixed<<setprecision(15)<<"Real Energy: "<<realenergy<<" Kcal/mol"<<"\n";
	end2 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds2 = end2 - start2;
    time_t end_time2 = std::chrono::system_clock::to_time_t(end2);
	cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds2.count() << " sec\n\n";

	/*Reciprocal Energy (k!=0) using the 2D FT and 1D FI method*/
	chrono::time_point<std::chrono::system_clock> start3, end3;
	start3 = chrono::system_clock::now();
	double recienergy_fft=reciprocal_fft(PosIons2, ion_charges, natoms, beta, boxcell, Kvec, grid ,order)*unitzer;
	cout<<fixed<<setprecision(15)<<"Reciprocal Energy FT: "<<recienergy_fft<<" Kcal/mol"<<"\n";
	end3 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds3 = end3 - start3;
    time_t end_time3 = std::chrono::system_clock::to_time_t(end3);
	cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds3.count() << " sec\n\n";

	/*Reciprocal Energy (k!=0) using the 3D FT, with the correction factor */
	chrono::time_point<std::chrono::system_clock> start8, end8;
	start8 = chrono::system_clock::now();
	double recienergy_pm=PM2DEwald(PosIons2, ion_charges, natoms, beta, boxcell, grid , Kvec, order)*unitzer;
	cout<<fixed<<setprecision(15)<<"Reciprocal Energy Correction: "<<recienergy_pm<<" Kcal/mol"<<"\n";
	end8 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds8 = end8 - start8;
    time_t end_time8 = std::chrono::system_clock::to_time_t(end8);
	cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds8.count() << " sec\n\n";

	/*Reciprocal Energy (k!=0) using the integral method*/
	chrono::time_point<std::chrono::system_clock> start4, end4;
	start4 = chrono::system_clock::now();
	double recienergy_ka=reciprocal_kawata(PosIons2, ion_charges, natoms, beta, boxcell, Kvec)*unitzer;
	cout<<fixed<<setprecision(15)<<"Reciprocal Energy Integral(k!=0): "<<recienergy_ka<<" Kcal/mol"<<"\n";
	end4 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds4 = end4- start4;
    time_t end_time4 = std::chrono::system_clock::to_time_t(end4);
	cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds4.count() << " sec\n\n";

	/*Reciprocal Energy (k==0)*/
	chrono::time_point<std::chrono::system_clock> start6, end6;
	start6 = chrono::system_clock::now();
	double recienergy_0=reci0(PosIons2, charge_prod, natoms, beta, boxcell)*unitzer;
	cout<<fixed<<setprecision(15)<<"Reciprocal Energy (k==0): "<<recienergy_0<<" Kcal/mol"<<"\n";
	end6 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds6 = end6- start6;
    time_t end_time6 = std::chrono::system_clock::to_time_t(end6);
	cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds6.count() << " sec\n\n";

	cout<<fixed<<setprecision(8)<< "Total: "<<recienergy_0+recienergy_fft<<" Kcal/mol"<<"\n";
	cout<<fixed<<setprecision(8)<< "Time: "<<elapsed_seconds6.count()+elapsed_seconds3.count()<<" sec"<<"\n";

	/*Reciprocal Energy (k==0) + Real Energy*/
	// chrono::time_point<std::chrono::system_clock> start7, end7;
	// start7 = chrono::system_clock::now();
	// vector<double> energy = realnreci0(PosIons2, charge_prod, natoms, beta, boxcell,cutoff);
	// cout<<fixed<<setprecision(15)<<"Reciprocal Energy (k==0): "<<energy[1]*unitzer<<" Kcal/mol"<<"\n";
	// cout<<fixed<<setprecision(15)<<"Real Energy : "<<energy[0]*unitzer<<" Kcal/mol"<<"\n";
	// end7 = chrono::system_clock::now();
	// chrono::duration<double> elapsed_seconds7 = end7- start7;
    // time_t end_time7 = std::chrono::system_clock::to_time_t(end7);
	// cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds7.count() << " sec\n\n";

	/*Testing the fiNUFFT library*/
	/*Without*/
	// chrono::time_point<std::chrono::system_clock> start7, end7;
	// start7 = chrono::system_clock::now();
	// double out = without(PosIons,natoms,boxcell,order);
	// end7 = chrono::system_clock::now();
	// chrono::duration<double> elapsed_seconds7 = end7- start7;
    // time_t end_time7 = std::chrono::system_clock::to_time_t(end7);
	// cout<<fixed<<setprecision(8)<< "Elapsed timer: " << elapsed_seconds7.count() << " sec\n\n";
	
	
	/*With*/
	// chrono::time_point<std::chrono::system_clock> start6, end6;
	// start6 = chrono::system_clock::now();
	
	// end6 = chrono::system_clock::now();
	// chrono::duration<double> elapsed_seconds6 = end6- start6;
    // time_t end_time6 = std::chrono::system_clock::to_time_t(end6);
	// cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds6.count() << " sec\n\n";

	// int mx=6,my=6;
	// double h = -0.0261105;
	// cout<<fixed<<setprecision(8)<<StructureFactor(mx,my,h,PosIons,ion_charges,natoms,boxcell,32,8)<<"\n";
	// cout<<fixed<<setprecision(8)<<StructureFactor2(mx,my,h,PosIons,ion_charges,natoms,boxcell,32,8)<<"\n";
	// cout<<fixed<<setprecision(8)<<StructureFactor3(mx,my,h,PosIons,ion_charges,natoms,boxcell,32,8)<<"\n";

	// delete dynamic variables 
	for(i=0;i<3;i++){
		delete [] boxcell[i]; 
	}

	for(i=0;i<natoms;i++){
		delete [] PosIons[i];
	}

	delete [] PosIons;
	delete [] PosIons2;
	delete [] boxcell;
	delete [] atomtype;
	delete [] natoms_type;
	delete [] ForceIons;
	delete [] vel;
	delete [] fixatoms;
    delete [] chg;
    delete [] ion_charges;
    delete [] charge_prod;

	return 0;
}
//end of main loop
