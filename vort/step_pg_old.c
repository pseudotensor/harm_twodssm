
#include "decs.h"

/* pressure and gravity step */
void step_pg()
{
	double Sa ;
	double XP,vxave,vyave ;
	void bound_var(REAL (* var)[NY+4], double t) ;
	int i,j ;

	/* first pressure gradient */
	if(fabs(g2 - 1.) > 1.e-6) {
		for(i=0;i<NX;i++)
		for(j=0;j<NY;j++) {
			Sa = 0.5*(S[i][j]+S[i-1][j]) ;
			vx[i][j] += -(dt/(dx*Sa))*(g2-1.)*(e[i][j]-e[i-1][j]) ;

			Sa = 0.5*(S[i][j]+S[i][j-1]) ;
			vy[i][j] += -(dt/(dy*Sa))*(g2-1.)*(e[i][j]-e[i][j-1]) ;
		}
	}
	else {
		for(i=0;i<NX;i++)
		for(j=0;j<NY;j++) {
			Sa = 0.5*(S[i][j]+S[i-1][j]) ;
			vx[i][j] += -(dt/(dx*Sa))*cs*cs*(S[i][j]-S[i-1][j]) ;

			Sa = 0.5*(S[i][j]+S[i][j-1]) ;
			vy[i][j] += -(dt/(dy*Sa))*cs*cs*(S[i][j]-S[i][j-1]) ;
		}
	}

	/* gravitational acceleration */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		vx[i][j] += -(dt/dx)*(pot[i][j]-pot[i-1][j]) ;
		vy[i][j] += -(dt/dy)*(pot[i][j]-pot[i][j-1]) ;
	}
	bound_var(vx,t) ;
	bound_var(vy,t) ;

	/* Coriolis force, tidal acceleration */
	for(i=0;i<NX;i++) {
#ifndef SWEEPYM
	XP = i*dx - 0.5*Lx ;
#endif
	for(j=0;j<NY;j++) {
		vyave = 0.25*(vy[i][j] + vy[i-1][j] +
			vy[i][j+1] + vy[i-1][j+1]) ;
		vx[i][j] += dt*2.*W*vyave ;
#ifndef SWEEPYM
		vx[i][j] += dt*2.*q*W*W*XP ;
#endif
	}
	}
	bound_var(vx,t) ;

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		vxave = 0.25*(vx[i][j] + vx[i+1][j] +
			vx[i][j-1] + vx[i+1][j-1]) ;
		vy[i][j] += dt*(-2.*W*vxave) ;
#ifdef SWEEPYM
		vy[i][j] += dt*q*W*vxave ;
#endif
	}
	bound_var(vy,t) ;

}

