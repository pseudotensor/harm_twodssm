/*
	produces an "r8" file.  
*/

#include "decs.h"

void image(FILE *fp,int im_var)
{
	REAL iq,liq,slope,intercept,lmax,lmin ;
	REAL Save,pv,h,hp ;
	double X ;
	int i,j,ilo ;
	
	if(im_var == 0) return ;

	lmax = -10. ;
	lmin = 10. ;

	if(im_var == 1) ilo = -4 ;
	else ilo = -2 ;	
	for(i=ilo;i<NX+2;i++) {
	X = (i+0.5)*dx - 0.5*Lx ;
        h = h0*(1. + eps*sin(kx0*X + phi)) ;
        hp = h0*eps*kx0*cos(kx0*X + phi) ;
	for(j=-2;j<NY+2;j++) {
		if(im_var == 1) {
			if(i==ilo) iq = 1. ;
			else iq = S[i][j] ;
		}
		if(im_var == 2) iq = log(e[i][j]) ;
		if(im_var == 3) {
			Save = 0.25*(S[i][j] + S[i-1][j] + S[i][j-1] + S[i-1][j-1]) ;
			pv = (    (2. - q)*W
			+ (vy[i][j] - vy[i-1][j])/dx
			- (vx[i][j] - vx[i][j-1])/dy
			)/Save ;
			iq = pv ;
		}
		if(im_var == 4) iq = vx[i][j] ;
		if(im_var == 5) iq = vy[i][j] ;
		if(iq > lmax) lmax = iq ;
		if(iq < lmin) lmin = iq ;
	}
	}

	slope = 256./(lmax - lmin) ;
	intercept = -slope*lmin ;

	for(j=NY+1;j>=-2;j--)
	for(i=ilo;i<NX+2;i++) {
	        X = (i+0.5)*dx - 0.5*Lx ;
	        h = h0*(1. + eps*sin(kx0*X + phi)) ;
	        hp = h0*eps*kx0*cos(kx0*X + phi) ;	
		if(im_var == 1) {
			if(i==ilo) iq = 1. ;
			else iq = S[i][j] ;
		}
		if(im_var == 2) iq = log(e[i][j]) ;
		if(im_var == 3) {
			Save = 0.25*(S[i][j] + S[i-1][j] + S[i][j-1] + S[i-1][j-1]) ;
			pv = (    (2. - q)*W
			+ (vy[i][j] - vy[i-1][j])/dx
			- (vx[i][j] - vx[i][j-1])/dy
			)/Save ;
			iq = pv ;
		}
		if(im_var == 4) iq = vx[i][j] ;
		if(im_var == 5) iq = vy[i][j] ;
		liq = slope*iq + intercept ;
		if(liq > 255.) liq = 255. ;
		if(liq < 0.) liq = 0. ;
		fprintf(fp,"%c",(char)((int)liq)) ;
	}
}
