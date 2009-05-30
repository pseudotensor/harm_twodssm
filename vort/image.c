/*
	produces an "r8" file.
*/

#include "decs.h"

void image(FILE *fp,int im_var)
{
	REAL iq,liq,slope,intercept,lmax,lmin ;
	REAL Save,pv,hmax,hmin,Smax,Smin,emax,emin,vymax,vymin,vxmax,vxmin ;
	REAL pvsq,pvav,sig ;
	double X ;
	int i,j,Ny=NY ;
	int ilo,nghost=NGHOST ;

	lmax = -10. ;
	lmin = 10. ;
	pvsq = 0. ;
	pvav = 0. ;
	for(j=NY-1+nghost;j>=-nghost;j--)
	for(i=ilo;i<NX+nghost;i++) {
		Save = 0.25*(S[i][j] + S[i-1][j] + S[i][j-1] + S[i-1][j-1]) ;
		pv = (    (2. - q)*W
			+ (vy[i][j] - vy[i-1][j])/dx
			- (vx[i][j] - vx[i][j-1])/dy
		)/Save ;
		iq = pv ;
		pvsq += pv*pv ;
		pvav += pv ;
		if(iq > lmax) lmax = iq ;
		if(iq < lmin) lmin = iq ;
	}
	pvsq /= NX*NY ;
	pvav /= NX*NY ;
	sig = sqrt( pvsq - pvav*pvav ) ;
	/*
	*/

	//fprintf(stderr,"max,min: %g %g\n",lmax,lmin) ;

	lmax = 0.5 + 3.*sig ;
	lmin = 0.5 - 3.*sig ;

	//fprintf(stderr,"lmax,lmin: %g %g\n",lmax,lmin) ;

	slope = 256./(lmax - lmin) ;
	intercept = -slope*lmin ;

	for(j=NY-1+nghost;j>=-nghost;j--)
	for(i=ilo;i<NX+nghost;i++) {
		Save = 0.25*(S[i][j] + S[i-1][j] + S[i][j-1] + S[i-1][j-1]) ;
		pv = (    (2. - q)*W
			+ (vy[i][j] - vy[i-1][j])/dx
			- (vx[i][j] - vx[i][j-1])/dy
		)/Save ;
		iq = pv ;

		liq = slope*iq + intercept ;
		if(liq > 255.) liq = 255. ;
		if(liq < 0.) liq = 0. ;
		fprintf(fp,"%c",(char)((int)liq)) ;
	}
}
