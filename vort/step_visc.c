
#include "decs.h"

void step_visc()
{
	static REAL (*visx)[NY+4],(*visy)[NY+4],(*sxx)[NY+4],(*syy)[NY+4],(*sxy)[NY+4],(*exx)[NY+4],
		(*eyy)[NY+4],(*exy)[NY+4],(*divv)[NY+4],(*diss)[NY+4] ;
	double dvxdx,dvydy,qlx,qvnrx,qly,qvnry ;
	double Sa,Sb,eb ;
	double H0,nu0 ;
	double sna,snb,dvxx,dvxy,dvyx,dvyy,phi,exyp ; 
	int i,j ;

	visx = work1 ;
	visy = work2 ;

	/* find viscous stresses */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		/* del v, at zone center */
		dvxdx = vx[i+1][j] - vx[i][j] ;
		dvydy = vy[i][j+1] - vy[i][j] ;
		
		/* von neumann,richtmyer viscosity */
		if(dvxdx  < 0) {
			qvnrx = nu_vnr*S[i][j]*dvxdx*dvxdx ;
		}
		else qvnrx = 0. ;
		if(dvydy  < 0) {
			qvnry = nu_vnr*S[i][j]*dvydy*dvydy ;
		}
		else qvnry = 0. ;

		visx[i][j] = qvnrx ;
		visy[i][j] = qvnry ;
	}
	bound_var(visx,t) ;
	bound_var(visy,t) ;

	/* update velocity, internal energy */
	if(fabs(g2 - 1.) > 1.e-6) {
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) {
			e[i][j] += 
				-dt*visx[i][j]*(vx[i+1][j] - vx[i][j])/dx 
				- dt*visy[i][j]*(vy[i][j+1] - vy[i][j])/dy ;
			if(e[i][j] < 0.) fprintf(stdout,"step_visc e < 0: %d %d %g %g\n",
					i,j,S[i][j],e[i][j]) ;
		}
	}

	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		Sa = 0.5*(S[i][j]+S[i-1][j]) ;
		vx[i][j] += -dt*(visx[i][j]-visx[i-1][j])/(dx*Sa) ;
		Sb = 0.5*(S[i][j]+S[i][j-1]) ;
		vy[i][j] += -dt*(visy[i][j]-visy[i][j-1])/(dy*Sb) ;
	}
	bound_var(vx,t) ;
	bound_var(vy,t) ;

}
