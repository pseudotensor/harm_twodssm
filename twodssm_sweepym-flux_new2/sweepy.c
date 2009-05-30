
#include "decs.h"

void sweepy()
{
	double transvarp,mdotp ;
	REAL (*dq)[NY+4] ;
	REAL (*mdot)[NY+4] ;
	REAL (*transvar)[NY+4] ;
	REAL (*px)[NY+4],(*py)[NY+4] ;
	int i,j,jp,jm ;

	px = work1 ;
	py = work2 ;
	mdot = work3 ;
	dq = work4 ;
	transvar = work5 ;

	/* transform to momenta */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		px[i][j] = 0.5*(S[i-1][j]+S[i][j])*vx[i][j] ;
		py[i][j] = 0.5*(S[i][j-1]+S[i][j])*vy[i][j] ;
	}
	/* no need to bound px,py */

	/* create dimensionless transport variable */
	for(i=-1;i<NX;i++)
	for(j=0;j<NY;j++) transvar[i][j] = vy[i][j]*dt/dy ;

	/* calculate mass flux */
	for(i=-1;i<NX;i++)
	for(j=0;j<NY;j++) {
		jp = j + 1 ;
		if(jp > NY-1) jp -= NY ;
		jm = j - 1 ;
		if(jm < 0)    jm += NY ;
		dq[i][j] = slope_lim(S[i][jm],S[i][j],S[i][jp]) ;
	}

	for(i=-1;i<NX;i++)
	for(j=0;j<NY;j++) {
		jm = j - 1 ;
		if(jm < 0)    jm += NY ;
		mdot[i][j] = (transvar[i][j] > 0.) ?
			transvar[i][j]*(S[i][jm] + (1. - transvar[i][j])*dq[i][jm]) :
			transvar[i][j]*(S[i][j] - (1. + transvar[i][j])*dq[i][j]) ;
	}

	/* get energy flux */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		jp = j + 1 ;
		if(jp > NY-1) jp -= NY ;
		jm = j - 1 ;
		if(jm < 0)    jm += NY ;
		dq[i][j] = slope_lim(e[i][jm]/S[i][jm],e[i][j]/S[i][j],e[i][jp]/S[i][jp]) ;
	}

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		jm = j - 1 ;
		if(jm < 0)    jm += NY ;
		fl[i][j] = (transvar[i][j] > 0) ?
			mdot[i][j]*(e[i][jm]/S[i][jm] + (1. - transvar[i][j])*dq[i][jm]) :
			mdot[i][j]*(e[i][j]/S[i][j] - (1. + transvar[i][j])*dq[i][j]) ;
	}

	/* step energy */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		jp = j + 1 ;
		if(jp > NY-1) jp -= NY ;
		e[i][j] += -(fl[i][jp]-fl[i][j]) ;
		if(e[i][j] < 0.) fprintf(stdout,"sweepy e < 0 %d %d %g %g\n",i,j,S[i][j],e[i][j]) ;
	}
	bound_var(e,t) ;

	/* step mass */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		jp = j + 1 ;
		if(jp > NY-1) jp -= NY ;
		S[i][j] += -(mdot[i][jp]-mdot[i][j]) ;
	}
	bound_var(S,t) ;

	/* get vx flux */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		jp = j + 1 ;
		if(jp > NY-1) jp -= NY ;
		jm = j - 1 ;
		if(jm < 0)    jm += NY ;
		dq[i][j] = slope_lim(vx[i][jm],vx[i][j],vx[i][jp]) ;
	}

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		jm = j - 1 ;
		if(jm < 0)    jm += NY ;
		transvarp = 0.5*(transvar[i][j]+transvar[i-1][j]) ;
		mdotp = 0.5*(mdot[i][j]+mdot[i-1][j]) ;
		fl[i][j] = (transvarp > 0) ?
			mdotp*(vx[i][jm] + (1. - transvarp)*dq[i][jm]) :
			mdotp*(vx[i][j] - (1. + transvarp)*dq[i][j]) ;
	}

	/* step vx */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		jp = j + 1 ;
		if(jp > NY-1) jp -= NY ;
		px[i][j] += -(fl[i][jp]-fl[i][j]) ;
	}

	/* get vy flux */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		jp = j + 1 ;
		if(jp > NY-1) jp -= NY ;
		jm = j - 1 ;
		if(jm < 0)    jm += NY ;
		dq[i][j] = slope_lim(vy[i][jm],vy[i][j],vy[i][jp]) ;
	}

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		jp = j + 1 ;
		if(jp > NY-1) jp -= NY ;
		transvarp = 0.5*(transvar[i][j]+transvar[i][jp]) ;
		mdotp = 0.5*(mdot[i][j]+mdot[i][jp]) ;
		fl[i][j] = (transvarp > 0.) ?
			mdotp*(vy[i][j] + (1. - transvarp)*dq[i][j]) :
			mdotp*(vy[i][jp] - (1. + transvarp)*dq[i][jp]) ;
	}

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		jm = j - 1 ;
		if(jm < 0)    jm += NY ;
		py[i][j] += -(fl[i][j]-fl[i][jm]) ;
	}

	/* return to velocities */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		vx[i][j] = 2.*px[i][j]/(S[i][j]+S[i-1][j]) ;
		vy[i][j] = 2.*py[i][j]/(S[i][j-1]+S[i][j]) ;
	}

	bound_var(vx,t) ;
	bound_var(vy,t) ;
}
