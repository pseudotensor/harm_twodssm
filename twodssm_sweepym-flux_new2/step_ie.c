
#include "decs.h"

void step_ie()
{
	double dv ;
	int i,j ;

	/* update internal energy */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		dv = (vx[i+1][j] - vx[i][j])/dx +
		     (vy[i][j+1] - vy[i][j])/dy ;

		e[i][j] *= (1. - 0.5*dt*(g2 - 1.)*dv)/(1. + 0.5*dt*(g2 - 1.)*dv) ;
		if(e[i][j] < 0.) {
			fprintf(stdout,"step_ie e < 0: %g %d %d %g %g\n",t,i,j,S[i][j],e[i][j]) ;
			//exit(0) ;
		}
	}
	bound_var(e,t) ;
}
