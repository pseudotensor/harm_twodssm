
#include "decs.h"

void sweepym()
{
	double del,X,XP ;
	int i,j ;

	for(i=0;i<NX;i++) {
		XP = i*dx - 0.5*Lx ;
		del = -q*W*XP*dt/dy ;
		remap_linear(vx,i,del) ;
	}

	for(i=0;i<NX;i++) {
                X = (i + .5)*dx - 0.5*Lx ;
		del = -q*W*X*dt/dy ;
		remap_linear(vy,i,del) ;
		remap_linear(e,i,del) ;
		remap_linear(S,i,del) ;
	}

	bound_var(vx,t+dt) ;
	bound_var(vy,t+dt) ;
	bound_var(e,t+dt) ;
	bound_var(S,t+dt) ;
}

