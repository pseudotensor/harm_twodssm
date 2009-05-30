
#include "decs.h"

void timestep() 
{
	int i,j ;
	double dt1,dt2,dt3,dt4,dt5,dt6,dt7 ;
	double idt1,idt2,idt3,idt4,idt5,idt6,idt7 ;
	REAL (*dtinv)[NY+4] ;
	double dtinv_max ;
	double dv,d ;
	double u,s,du,ds ;
	static double dtlast ;
	double mymax(double a, double b) ;
	void dump() ;

	dtinv_max = 0. ;
	dtlast = dt ;

	if(dx > dy) d = dy ;
	else d = dx ;

	/* SSBC timestep criterion */

	dtinv = work1 ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {

		if(S[i][j] < 0) {
//			dump(stdout) ;
			fprintf(stdout,"negative density error\n") ;
                        fprintf(stdout,"%g %d %d %10.5g %10.5g\n",t,i,j,S[i][j],e[i][j]) ;
			exit(4) ;
		}

		/* x-velocity */
		idt2 = (fabs(vx[i][j]) + cs)/dx ;

		/* y-velocity */
		idt3 = (fabs(vy[i][j]) + cs)/dy ;

		dtinv[i][j] = mymax(idt2,idt3) ;

	}
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		if(dtinv[i][j] > dtinv_max) {
			dtinv_max = dtinv[i][j] ;
			i_crit_dt = i ;
			j_crit_dt = j ;
		}
	}


	dt = cour/dtinv_max ;

	//fprintf(stderr,"%g %g %g\n",cs,dtinv_max,dt) ;

	/* don't increase timestep by too much */
	if(dt > 1.3*dtlast) dt = 1.3*dtlast ;

	/* don't step beyond next dumping time */
	if(dt > 0.5*DTd) dt = 0.5*DTd ;

	/* don't step beyond end of run */
	if(t + dt > tf) dt = tf - t ;
}

double mymax(double a, double b)
{
	if(a > b) return(a) ;
	else return(b) ;
}
