
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define REAL double
#define NX 128
#define NY 128

extern REAL S_matrix[NX+5][NY+4] ;
extern REAL Stmp_matrix[NX+5][NY+4] ;
extern REAL e_matrix[NX+4][NY+4] ;
extern REAL vx_matrix[NX+4][NY+4] ;
extern REAL vy_matrix[NX+4][NY+4] ;
extern REAL pot_matrix[NX+4][NY+4] ;
extern REAL fl_matrix[NX+4][NY+4] ;
extern REAL N_matrix[NX+4][NY+4] ;
extern REAL work1_matrix[NX+4][NY+4] ;
extern REAL work2_matrix[NX+4][NY+4] ;
extern REAL work3_matrix[NX+4][NY+4] ;
extern REAL work4_matrix[NX+4][NY+4] ;
extern REAL work5_matrix[NX+4][NY+4] ;
extern REAL work6_matrix[NX+4][NY+4] ;

extern REAL (* S)[NY+4] ;
extern REAL (* Stmp)[NY+4] ;
extern REAL (* e)[NY+4] ;
extern REAL (* vx)[NY+4] ;
extern REAL (* vy)[NY+4] ;
extern REAL (* pot)[NY+4] ;
extern REAL (* fl)[NY+4] ;
extern REAL (* N)[NY+4] ;
extern REAL (* work1)[NY+4] ;
extern REAL (* work2)[NY+4] ;
extern REAL (* work3)[NY+4] ;
extern REAL (* work4)[NY+4] ;
extern REAL (* work5)[NY+4] ;
extern REAL (* work6)[NY+4] ;

extern double cour ;
extern double cs ;
extern double g2,deltag ;
extern double dx,dy ;
extern double dt ;
extern double t,tf ;
extern double nu_vnr ;
extern double nu_l ;
extern double G ;
extern double cool_fac ;
extern double mu ;

extern int i_crit_dt ;
extern int j_crit_dt ;
extern int i_crit_tc ;
extern int j_crit_tc ;
extern int dt_crit ;

extern double DTd ;
extern double DTi ;
extern double DTl ;
extern double res ;
extern double nu_sh ;

extern double W ;
extern double omega ;
extern double q ;
extern double Lx ;
extern double Ly ;
extern double K0 ;
extern double alphah,alphag,cs0,u0,S0 ;
extern double kx0,ky0 ;
extern double kmax,ksqmax ;
extern double fx,fy ;
extern double Qeq ;
extern double pamp ;
extern int is_restart ;
extern double h0,phi,eps ;
extern double N2min,Xmin ;
extern double dampinit,damp ;
extern double a,ar ;

extern int dump_cnt ;
extern int im_cnt ;

#define SMALL	1.e-20
#define SMIN	1.e-6	/* minimum density */
#define EMIN	1.e-6	/* minimum internal energy */

#define CONSTVISC	1
#define VARVISC		0

/* subroutines */

double pow(double x, double y) ;

int my_nint(double dum) ;

float *vector(int xl, int xh) ;

float **matrix(int xl, int xh, int yl, int yh) ;

double **dmatrix(int xl, int xh, int yl, int yh) ;

double ranc(int seed) ;

void bound_var(REAL (* var)[NY+4], double tcurr) ;
void consistent_transport(REAL (*var)[NY+4], int i, REAL *massvar, REAL *specvar, REAL *fluxvar, int fr) ;
void diag(int call_code) ;
void dump(FILE *fp) ;
void fft(int n,double *real,double *imag) ;
void fft_grid(REAL **inr,REAL **ini,REAL **outr,REAL **outi,int nx,int ny,int dir) ;
void flux_calc(REAL *var, REAL *trnsvar, REAL *mdotvar, REAL *flux, int offset) ;
void flux_update(REAL *fluxvar, REAL *flux, int offset) ;
void image(FILE *fp,int im_var) ;
void iremap(REAL *var, int shift) ;
void dperturb(double amp) ;
void potsolv(REAL (* dens)[NY+4], REAL (* psi)[NY+4]) ;
void remap_flux(REAL (*var)[NY+4], int i, int shift, double del) ;
void remap_linear(REAL (*var)[NY+4], int i, double del) ;
void remap_fourier(REAL (*var)[NY+4], int i, double del) ;
void row_unroll(REAL *var, REAL del) ;
double slope_lim(double y1,double y2,double y3) ;
void tfft(REAL *inr1,REAL *ini1,REAL *outr1,REAL *outi1, int n, int dir) ;
void unroll(REAL (*var)[NY+4], REAL delt) ;



