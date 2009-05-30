
/* function for shifting a column of data forward by some amount */

#include "decs.h"

void remap_flux(REAL (*var)[NY+4], int i, int shift, double del)
{
	int j,jp,jm,idel,offset ;
	double fdel ;
	REAL dq[NY],flux[NY],mdot[NY],trnsvar[NY] ;
	REAL massvar[NY],specvar[NY],fluxvar[NY] ;

        /* separate into integral & fraction part of shift */
        idel = my_nint(del) ;
        fdel = del - (double)idel ;

	offset = (var==vy) ? 1 : 0 ;

	consistent_transport(var,i+shift,massvar,specvar,fluxvar,1) ;

	/* calculate mass flux */
	for(j=0;j<NY;j++) trnsvar[j] = fdel ;
	flux_calc(massvar,trnsvar,trnsvar,mdot,offset) ;

	/* calculate flux */
	if(var==S||var==Stmp) for(j=0;j<NY;j++) flux[j] = mdot[j] ;
	else flux_calc(specvar,trnsvar,mdot,flux,offset) ;

	/* update transport variable */
	flux_update(fluxvar,flux,offset) ;

	/* do integral part of shift */
        iremap(fluxvar,idel) ;

	/*if(var==vx || var==vy) {
		flux_update(massvar,mdot,offset) ;
		iremap(massvar,idel) ;
	}*/

	consistent_transport(var,i,massvar,specvar,fluxvar,-1) ;
}


void iremap(REAL *var, int idel)
{
        static REAL tmp[NY] ;
	int k,kp ;

	for(k=0;k<NY;k++) tmp[k] = var[k] ;

        for(k=0;k<NY;k++) {
                kp = k + idel ;
                while(kp < 0) kp += NY ;
                while(kp > NY-1) kp -= NY ;
                var[kp] = tmp[k] ;
	}
}


void consistent_transport(REAL (*var)[NY+4], int i, REAL *massvar, REAL *specvar, REAL *fluxvar, int fr)
{
	int j,jm ;

	if(fr==1) {
		if(var==vx) for(j=0;j<NY;j++) {
			massvar[j] = 0.5*(S[i-1][j]+S[i][j]) ;
			specvar[j] = var[i][j] ;
			fluxvar[j] = massvar[j]*var[i][j] ;
		}
		else if(var==vy) for(j=0;j<NY;j++) {
			massvar[j] = 0.5*(S[i][j-1]+S[i][j]) ;
			specvar[j] = var[i][j] ;
			fluxvar[j] = massvar[j]*var[i][j] ;
		}
		else if(var==S||var==Stmp) for(j=0;j<NY;j++) {
			massvar[j] = fluxvar[j] = S[i][j] ;
		}
		else for(j=0;j<NY;j++) {
			massvar[j] = S[i][j] ;
			specvar[j] = var[i][j]/massvar[j] ;
			fluxvar[j] = var[i][j] ;
		}
	}

	else if(fr==-1) {
		if(var==vx) for(j=0;j<NY;j++) var[i][j] = 2.0*fluxvar[j]/(Stmp[i-1][j]+Stmp[i][j]) ;
		else if(var==vy) for(j=0;j<NY;j++) {
			jm = j - 1 ;
			if(jm < 0) jm +=NY ;
			var[i][j] = 2.0*fluxvar[j]/(Stmp[i][jm]+Stmp[i][j]) ;
		}
		//if(var==vx || var==vy) for(j=0;j<NY;j++) var[i][j] = fluxvar[j]/massvar[j] ;
  		else for(j=0;j<NY;j++) var[i][j] = fluxvar[j] ;
	}
	else exit(0) ;
}


void flux_calc(REAL *var, REAL *trnsvar, REAL *mdotvar, REAL *flux, int offset)
{
	int k,km,kp ;
	REAL dq[NY] ;

	for(k=0;k<NY;k++) {
		kp = k + 1 ;
		if(kp > NY-1) kp -= NY ;
		km = k - 1 ;
		if(km < 0)    km += NY ;
		dq[k] = slope_lim(var[km],var[k],var[kp]) ;
	}

	for(k=0;k<NY;k++) {
		kp = k + offset ;
		if(kp > NY-1) kp -= NY ;
		km = kp - 1 ;
		if(km < 0)    km += NY ;
		flux[k] = (trnsvar[k] > 0.) ?
			mdotvar[k]*(var[km] + (1. - trnsvar[k])*dq[km]) :
			mdotvar[k]*(var[kp] - (1. + trnsvar[k])*dq[kp]) ;
	}
}


void flux_update(REAL *fluxvar, REAL *flux, int offset)
{
	int k,km,kp ;

	for(k=0;k<NY;k++) {
		km = k - offset ;
		if(km < 0)    km += NY ;
		kp = km + 1 ;
		if(kp > NY-1) kp -= NY ;
		fluxvar[k] += -(flux[kp]-flux[km]) ;
	}
}


double slope_lim(double y1,double y2,double y3)
{
        double Dqm,Dqp,s ;

        /* van leer slope limiter */
	/* N.B.: factor of 2 is missing because factor is
	   absorbed into expressions for fluxes above */
	Dqp = (y3 - y2) ;
	Dqm = (y2 - y1) ;
	s = Dqm*Dqp ;
	if(s <= 0.) return 0. ;
	else return(s/(Dqm+Dqp)) ;

}
