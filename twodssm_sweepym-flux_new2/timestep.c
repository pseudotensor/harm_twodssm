
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
	void dump() ;

	dtinv_max = 0. ;
	dtlast = dt ;

	if(dx > dy) d = dy ;
	else d = dx ;

	/* SSBC timestep criterion */
	/*
	dt4 = dy/(q*W*Lx + SMALL) ;
	idt4 = q*W*Lx/dy ;
	*/
	idt4 = 0. ;

	dtinv = work1 ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		u = e[i][j] ;
		if(u < 0) {
//			dump(stdout) ;
			fprintf(stdout,"negative internal energy error\n") ;
			fprintf(stdout,"%g %d %d %10.5g %10.5g\n",t,i,j,S[i][j],e[i][j]) ;
			exit(3) ;
		}
		if(S[i][j] < 0) {
//			dump(stdout) ;
			fprintf(stdout,"negative density error\n") ;
                        fprintf(stdout,"%g %d %d %10.5g %10.5g\n",t,i,j,S[i][j],e[i][j]) ;
			exit(4) ;
		}

		/* sound speed */
		if(fabs(g2 - 1.) > SMALL) cs = sqrt(g2*(g2-1.)*u/S[i][j]) ;
		idt1 = cs/d ;
		work2[i][j] = 1 ;

		/* x-velocity */
		idt2 = vx[i][j]/dx ;
		if(idt2 > idt1) work2[i][j] = 2 ; 

		/* y-velocity */
		idt3 = vy[i][j]/dy ;
		if(idt3 > idt1 && idt3 > idt2) work2[i][j] = 3 ;

		/* viscosities */
		idt5 = 4.*nu_l*cs/d ;
		if(idt5 > idt1 && idt5 > idt2 && idt5 > idt3) work2[i][j] = 5 ;

		dv = vx[i+1][j]-vx[i][j] ;
		idt6 = 4.*nu_vnr*dv/dx ;
		if(idt6 > idt1 && idt6 > idt2
 			&& idt6 > idt3 && idt6 > idt5) work2[i][j] = 6 ;

		dv = vy[i][j+1]-vy[i][j] ;
		idt7 = 4.*nu_vnr*dv/dy ;
		if(idt7 > idt1 && idt7 > idt2 && idt7 > idt3
                        && idt7 > idt5 && idt7 > idt6) work2[i][j] = 7 ;

		dtinv[i][j] = 
			idt1*idt1 +
			idt2*idt2 +
			idt3*idt3 +
			idt4*idt4 +
			idt5*idt5 +
			idt6*idt6 +
			idt7*idt7 ;

		/*
		fprintf(stderr,"%g %g %g %g %g %g %g\n",idt1,idt2,idt3,idt4,
			idt5,idt6,idt7) ;
		*/
	}
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		if(dtinv[i][j] > dtinv_max) {
			dtinv_max = dtinv[i][j] ;
			i_crit_dt = i ;
			j_crit_dt = j ;
			dt_crit = work2[i][j] ;
		}
	}

	dt = cour/sqrt(dtinv_max) ;

	/* don't increase timestep by too much */
	if(dt > 1.3*dtlast) dt = 1.3*dtlast ;

	/* don't step beyond next dumping time */
	if(dt > 0.5*DTd) dt = 0.5*DTd ;

	/* don't step beyond end of run */
	if(t + dt > tf) dt = tf - t ;
}
