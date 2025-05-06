///declare all functions here

void print_lammps_input_file(double *PosIons, double *chg, int natoms, double **boxcell, int n_atomtype, int *natoms_type, string *atomtype, int printtrj, int MDstep, char printmode, string filename);

void print_coor(double *PosIons, int natoms, double **boxcell, int n_atomtype, int *natoms_type, string *atomtype, int printtrj, int MDstep,char printmode, string filename);

void print_carcoor(double *PosIons, int natoms, double **boxcell, int n_atomtype, int *natoms_type, string *atomtype, int printtrj, int MDstep,char printmode, string filename);

void printCoor(double *PosIons, int natoms, string *type);

void printFor(double *ForceIons, int natoms, string *type);

void printVel(double *Vel, int natoms, string *type);

double self(int n_atomtype, int *natoms_type, double *chargs, double betaa);

double real(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, double cutoff);

double realExact(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, double cutoff);

double F_0(double val);

double F_kl(double *ri, double *rj, double sigma, double psi, double beta, bool same_r, double **box);

double F_kl_I(double delX, double delY, double delZ, double sigma, double psi, double norm_sigma_psi, double beta, double **box);

double F_kl_0(double norm_sigma_psi, double beta);

double reciprocal_n2(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, int *K);

double reciprocal_modified(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, int *K);

double reciprocal_kawata(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, int* K);

double integrand_reciprocal(double h, void *params);

double reci0(double *PosIons2, double *ion_charges, int natoms, double betaa, double **box);

vector<double> realnreci0(double *PosIons2, double *charge_prod, int natoms, double betaa, double **box, double cutoff);

double dist(double *PosIons2, int atom1, int atom2, double **box);

template<typename T>
double dotProduct(T v1, T v2, size_t n) ;

template<typename T1>
void crossProduct(T1 v_A, T1 v_B, double *out);

double M_n(double u, int n);

complex<double>Coeff(double v, double w);

template<typename T2>
double RelError(T2 a, T2 b);

template<typename T3>
double percentReduction(T3 newValue, T3 oldValue);

double* linspace(int a, int b, int num);

double reciprocal_ft_integrand(double h, void *params);

double reciprocal_fft(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, int* K, int *Grid, int *n);

double PM2DEwald(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, int* Grid, int *M, int* n);

double tophat(double s, double gamma, double L);

double C_integrand(double s, void *params);

double C(double t, double n, double gamma, double L);

double constant_integrand(double t, void *params);

double constantterm(int kx, int ky, int kz, double lx, double ly, double lz, double beta, double gamma);