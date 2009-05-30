
#include "decs.h"

void sweepx()
{
	double transvarp,mdotp ;
	REAL (*dq)[NY+4] ;
	REAL (*mdot)[NY+4] ;
	REAL (*transvar)[NY+4] ;
	REAL (*px)[NY+4],(*py)[NY+4] ;
	int i,j,jm ;

	/* set auxiliary array names */
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
	} /* ... no need to bound px,py */

	/* create dimensionless transport variable */
	for(i=-2;i<NX+2;i++)
	for(j=0;j<NY;j++) transvar[i][j] = vx[i][j]*dt/dx ;
	//bound_var(vx,t) ;

	/* calculate mass flux */
	for(i=-2;i<NX+1;i++)
	for(j=0;j<NY;j++) dq[i][j] = slope_lim(S[i-1][j],S[i][j],S[i+1][j]) ;

	for(i=-1;i<NX+1;i++)
	for(j=0;j<NY;j++) {
		mdot[i][j] = (transvar[i][j] > 0.) ?
			transvar[i][j]*(S[i-1][j] + (1. - transvar[i][j])*dq[i-1][j]) :
			transvar[i][j]*(S[i][j] - (1. + transvar[i][j])*dq[i][j]) ;
	}

	/* get energy flux */
	for(i=-1;i<NX+1;i++)
	for(j=0;j<NY;j++) dq[i][j] = slope_lim(e[i-1][j]/S[i-1][j],
			e[i][j]/S[i][j],e[i+1][j]/S[i+1][j]) ;

	for(i=0;i<NX+1;i++)
	for(j=0;j<NY;j++) {
		fl[i][j] = (transvar[i][j] > 0.) ?
			mdot[i][j]*(e[i-1][j]/S[i-1][j] + (1. - transvar[i][j])*dq[i-1][j]) :
			mdot[i][j]*(e[i][j]/S[i][j] - (1. + transvar[i][j])*dq[i][j]) ;
	}

	/* step energy */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		e[i][j] += -(fl[i+1][j]-fl[i][j]) ;
		if(e[i][j] < 0.) fprintf(stdout,"sweepx e < 0 %d %d %g %g\n",i,j,S[i][j],e[i][j]) ;
	}
	bound_var(e,t) ;

	/* step mass */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		S[i][j] += -(mdot[i+1][j]-mdot[i][j]) ;
	}
	bound_var(S,t) ;

	/* get vx flux */
	for(i=-1;i<NX+1;i++)
	for(j=0;j<NY;j++) dq[i][j] = slope_lim(vx[i-1][j],vx[i][j],vx[i+1][j]) ;

	for(i=-1;i<NX;i++)
	for(j=0;j<NY;j++) {
		transvarp = 0.5*(transvar[i+1][j]+transvar[i][j]) ;
		mdotp = 0.5*(mdot[i+1][j]+mdot[i][j]) ;
		fl[i][j] = (transvarp > 0.) ?
			mdotp*(vx[i][j] + (1. - transvarp)*dq[i][j]) :
			mdotp*(vx[i+1][j] - (1. + transvarp)*dq[i+1][j]) ;
	}

	/* step vx */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) px[i][j] += -(fl[i][j]-fl[i-1][j]) ;

	/* get vy flux */
	for(i=-1;i<NX+1;i++)
	for(j=0;j<NY;j++) dq[i][j] = slope_lim(vy[i-1][j],vy[i][j],vy[i+1][j]) ;

	/* step vy */
	for(i=0;i<NX+1;i++)
	for(j=0;j<NY;j++) {
		jm = j - 1 ;
		if(jm < 0)    jm += NY ;
		mdotp = 0.5*(mdot[i][j]+mdot[i][jm]) ;
		transvarp = 0.5*(transvar[i][j]+transvar[i][jm]) ;
		fl[i][j] = (transvarp > 0.) ?
			mdotp*(vy[i-1][j] + (1. - transvarp)*dq[i-1][j]) :
			mdotp*(vy[i][j] - (1. + transvarp)*dq[i][j]) ;
	}

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) py[i][j] += -(fl[i+1][j]-fl[i][j]) ;

	/* return to velocities */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		vx[i][j] = 2.*px[i][j]/(S[i][j]+S[i-1][j]) ;
		vy[i][j] = 2.*py[i][j]/(S[i][j]+S[i][j-1]) ;
	}
	bound_var(vx,t) ;
	bound_var(vy,t) ;
}