
#include "decs.h"
#define KMAX	32
#define AMPF	0.2

void init_spec()
{
	void gauss_gen() ;
	double norm,**dmatrix();
	double **az ;
	double n ;
	int i,j,ip,jp ;

	az = dmatrix(0,NX-1,0,NY-1) ;

	/* generate random fields */
	n = 0. ;	/* E_k \sim k^(-n) */
	gauss_gen(az,-0.5*(n + 3.)) ;

	/* difference potential to get divergence-free velocity */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		ip = i+1 ;
		jp = j+1 ;
		if(ip == NX) ip = 0 ;
		if(jp == NY) jp = 0 ;
		vx[i][j] = (az[i][jp] - az[i][j])/dy ;
		vy[i][j] = -(az[ip][j] - az[i][j])/dx ;
	}

	/* normalize */
	norm = 0. ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		/* assumes density = 1 */
		norm += 0.5*(
			vx[i][j]*vx[i][j] +
			vy[i][j]*vy[i][j]) ;
	}

	norm = sqrt(fabs(AMPF)/(norm*dx*dy)) ;

	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		vx[i][j] *= norm ;
		vy[i][j] *= norm ;
	}

	free_dmatrix(az,0,NX-1,0,NY-1) ;
}

/* add Gaussian random field to matrix f, with power spectrum
   amplitude norm */
void gauss_gen(f,sl)
double **f ;
double sl ; 
{
	float **dfr,**dfi ;
	float **matrix() ;
	double kx,ky,phase,amp,powsp(),ranc() ;
	void fft_grid() ;
	int i,j,ip,jp ;

	dfr = matrix(0,NX-1,0,NY-1) ;
	dfi = matrix(0,NX-1,0,NY-1) ;

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) dfr[i][j] = dfi[i][j] = 0. ;

	/* generate complex random field; use only real part */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		if(i > NX/2) ip = NX-i ;
		else ip = i ;
		if(j > NY/2) jp = NY-j ;
		else jp = j ;

		if((ip != 0 || jp != 0) && ip*ip+jp*jp < KMAX*KMAX) {
			phase = 2.*M_PI*ranc(0) ;
			amp = powsp((double)ip,(double)jp,sl) ;
		}
		else {
			phase = 0. ;
			amp = 0. ;
		}

		dfr[i][j] = (float)amp*cos(phase) ;
		dfi[i][j] = (float)amp*sin(phase) ;
	}

	/* fft */
	fft_grid(dfr,dfi,dfr,dfi,NX,NY,-1) ;

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) f[i][j] = (double)dfr[i][j] ;

	free_matrix(dfr,0,NX-1,0,NY-1) ;
	free_matrix(dfi,0,NX-1,0,NY-1) ;

	return ;
}

double powsp(kx,ky,sl)
double kx,ky,sl ;
{
	double k2,p ;
	double ranc(),x1,x2 ;


	k2 = kx*kx + ky*ky ;
	p = pow(k2,0.5*sl) ;

	x1 = ranc(0) ;
	x2 = ranc(0) ;	/* wasteful */
	p *= sqrt(-2.*log(x1))*cos(2.*M_PI*x2) ;


	return(p) ;
}
