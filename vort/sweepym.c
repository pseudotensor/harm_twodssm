
#include "decs.h"

void sweepym()
{
	double del,X,XP ;
	int i,j ;
	
#ifdef FLUX
	for(i=-1;i<NX;i++) {
		XP = i*dx - 0.5*Lx ;
		del = -q*W*XP*dt/dy ;
		remap_flux(Stmp,i,0,del) ;
	}
#endif

	for(i=0;i<NX;i++) {
		XP = i*dx - 0.5*Lx ;
		del = -q*W*XP*dt/dy ;
#if defined(LINEAR)
		remap_linear(vx,i,del) ;
#elif defined(FOURIER)
		remap_fourier(vx,i,del) ;
#elif defined(FLUX)
		remap_flux(vx,i,0,del) ;
#endif
	}

#ifdef FLUX
	for(i=0;i<NX;i++) {
		X = (i + .5)*dx - 0.5*Lx ;
		del = -q*W*X*dt/dy ;
		remap_flux(Stmp,i,0,del) ;
	}
#endif

	for(i=0;i<NX;i++) {
                X = (i + .5)*dx - 0.5*Lx ;
		del = -q*W*X*dt/dy ;
#if defined(LINEAR)
		remap_linear(vy,i,del) ;
		remap_linear(e,i,del) ;
		remap_linear(S,i,del) ;
#elif defined(FOURIER)
		remap_fourier(vy,i,del) ;
		remap_fourier(e,i,del) ;
		remap_fourier(S,i,del) ;
#elif defined(FLUX)
		remap_flux(vy,i,0,del) ;
		remap_flux(e,i,0,del) ;
		remap_flux(S,i,0,del) ;
#endif
	}

	bound_var(vx,t+dt) ;
	bound_var(vy,t+dt) ;
	bound_var(e,t+dt) ;
	bound_var(S,t+dt) ;
}

