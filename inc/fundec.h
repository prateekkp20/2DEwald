///declare all functions here

void print_lammps_input_file(double **PosIons, float *chg, int natoms, double **boxcell, int n_atomtype, int *natoms_type, string *atomtype, int printtrj, int MDstep, char printmode, string filename);

void print_coor(double **PosIons, int natoms, double **boxcell, int n_atomtype, int *natoms_type, string *atomtype, int printtrj, int MDstep,char printmode, string filename);

void print_carcoor(double **PosIons, int natoms, double **boxcell, int n_atomtype, int *natoms_type, string *atomtype, int printtrj, int MDstep,char printmode, string filename);

void printCoor(double **PosIons, int natoms, string *type);

void printFor(double **ForceIons, int natoms, string *type);

void printVel(double **Vel, int natoms, string *type);

void printprobVel(double **vel, int natoms, double *mass, float Temp);

double self(int n_atomtype, int *natoms_type, float *chargs, double betaa);

double real(double *PosIons2, double *charge_prod, int natoms, double betaa, double **box, double cutoff);

double F_0(double val);

double F_0_New(double DelZ, double DelZ2, double beta);

double F_kl(double *ri, double *rj, double sigma, double psi, double beta, bool same_r, double **box);

double reciprocal_n2(double **PosIons, double *charge_prod,  double *ion_charges, int natoms, double betaa, double **box, int K);

double reciprocal_kawata(double **PosIons, double *ion_charges, int natoms, double betaa, double **box, int K);

double integrand_reciprocal(double h, void *params);

double reci0(double *PosIons2, double *charge_prod, int natoms, double betaa, double **box);

vector<double> realnreci0(double *PosIons2, double *charge_prod, int natoms, double betaa, double **box, double cutoff);

double dist(double *PosIons2, int atom1, int atom2, double **box);

template<typename T>
double dotProduct(T v1, T v2, size_t n) ;

template<typename T1>
void crossProduct(T1 v_A, T1 v_B, double *out);

double M_n(double u, int n);

complex<double>Coeff(double v, double w);

template<typename T2>
double error(T2 a, T2 b);

template<typename T3>
double percentReduction(T3 newValue, T3 oldValue);

template<typename T4>
T4* linspace(T4 a, T4 b, T4 num);

double reciprocal_fastft_integrand(double h, void *params);

double reciprocal_ft_integrand(double h, void *params);

double reciprocal_fft(double **PosIons, double *ion_charges, int natoms, double betaa, double **box, int K, int Grid, int n);

double reciprocal_pppm(double **PosIons, double *ion_charges, int natoms, double betaa, double **box, int K, int Grid[], int n[]);

double reciprocal_pppm_chebyshev(double **PosIons, double *ion_charges, int natoms, double betaa, double **box, int K, int Grid[], int n[]);