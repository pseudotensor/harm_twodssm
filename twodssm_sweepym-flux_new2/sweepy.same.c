
#include "decs.h"

void sweepy()
{
	double trnsvarp,mdotp ;
	REAL (*dq)[NY+4] ;
	REAL (*mdotm)[NY+4] ;
	REAL (*trnsvarm)[NY+4] ;
	REAL (*px)[NY+4],(*py)[NY+4] ;
	REAL trnsvar[NY],trnsvarpx[NY],trnsvarpy[NY] ;
	REAL mdot[NY],mdotpx[NY],mdotpy[NY] ;
	REAL massvar[NY],specvar[NY],fluxvar[NY],flux[NY] ;
	int i,j,jp,jm ;

	px = work1 ;
	py = work2 ;
	mdotm = work3 ;
	dq = work4 ;
	trnsvarm = work5 ;

	/* create dimensionless transport variable */
	for(i=-1;i<NX;i++)
	for(j=0;j<NY;j++) {
		trnsvarm[i][j] = vy[i][j]*dt/dy ;
		Stmp[i][j] = S[i][j] ;
	}

	/* calculate mass flux */
	for(i=0;i<NX;i++) {
		flux_calc(S[i],trnsvarm[i],trnsvarm[i],mdot,0) ;
		for(j=0;j<NY;j++) mdotm[i][j] = mdot[j] ;
	}

	for(i=0;i<NX;i++) {
		consistent_transport(e,i,massvar,specvar,fluxvar,1) ;
		//flux_calc(massvar,trnsvarm[i],trnsvarm[i],mdot,0) ;
		flux_calc(specvar,trnsvarm[i],mdotm[i],flux,0) ;
		flux_update(fluxvar,flux,0) ;
		consistent_transport(e,i,massvar,specvar,fluxvar,-1) ;
		consistent_transport(vx,i,massvar,specvar,fluxvar,1) ;
		for(j=0;j<NY;j++) {
			trnsvarpx[j] = 0.5*(trnsvarm[i][j]+trnsvarm[i-1][j]) ;
			mdotpx[j] = 0.5*(mdotm[i][j]+mdotm[i-1][j]) ;
		}
		//flux_calc(massvar,trnsvarpx,trnsvarpx,mdotpx,0) ;
		flux_calc(vx[i],trnsvarpx,mdotpx,flux,0) ;
		flux_update(fluxvar,flux,0) ;
		for(j=0;j<NY;j++) px[i][j] = fluxvar[j] ;
		consistent_transport(vy,i,massvar,specvar,fluxvar,1) ;
		for(j=0;j<NY;j++) {
			jp = j + 1 ;
			if(jp > NY-1) jp -= NY ;
			trnsvarpy[j] = 0.5*(trnsvarm[i][j]+trnsvarm[i][jp]) ;
			mdotpy[j] = 0.5*(mdotm[i][j]+mdotm[i][jp]) ;
		}
		//flux_calc(massvar,trnsvarpy,trnsvarpy,mdotpy,1) ;
		flux_calc(vy[i],trnsvarpy,mdotpy,flux,1) ;
		flux_update(fluxvar,flux,1) ;
		for(j=0;j<NY;j++) py[i][j] = fluxvar[j] ;
	}

	for(i=0;i<NX;i++) flux_update(S[i],mdotm[i],0) ;

	bound_var(e,t) ;
	bound_var(S,t) ;

	/* return to velocities */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		vx[i][j] = 2.*px[i][j]/(S[i][j]+S[i-1][j]) ;
		vy[i][j] = 2.*py[i][j]/(S[i][j-1]+S[i][j]) ;
	}

	bound_var(vx,t) ;
	bound_var(vy,t) ;
}
