
/* enforce periodic boundary conditions */

#include "decs.h"

/* enforce boundary conditions for variable var at time tcurr */
void bound_var(REAL (* var)[NY+4], double tcurr)
{
	double del ;
	int i,j,idel,offset ;
	REAL massvar[NY],specvar[NY],fluxvar[NY] ;

	/* x-boundary */
	del = -q*W*Lx*tcurr/dy ;

	for(i=NX;i<=NX+1;i++) {
	if(q==0.) for(j=0;j<NY;j++) var[i][j] = var[i-NX][j] ;
	else {
#if defined(LINEAR)
		for(j=0;j<NY;j++) var[i][j] = var[i-NX][j] ;
		remap_linear(var,i,del) ;
#elif defined(FOURIER)
		for(j=0;j<NY;j++) var[i][j] = var[i-NX][j] ;
		remap_fourier(var,i,del) ;
#elif defined(FLUX)
		if(var==vx||var==vy) remap_flux(Stmp,i,-NX,del) ;
		if(var==vx) remap_flux(Stmp,i-1,-NX,del) ;
		remap_flux(var,i,-NX,del) ;
#endif
	}
	}
	for(i=-2;i<=-1;i++) {
	if(q==0.) for(j=0;j<NY;j++) var[i][j] = var[i+NX][j] ;
	else {
#if defined(LINEAR)
		for(j=0;j<NY;j++) var[i][j] = var[i+NX][j] ;
		remap_linear(var,i,-del) ;
#elif defined(FOURIER)
		for(j=0;j<NY;j++) var[i][j] = var[i+NX][j] ;
		remap_fourier(var,i,-del) ;
#elif defined(FLUX)
		if(var==vx||var==vy) remap_flux(Stmp,i,NX,-del) ;
		if(var==vx) remap_flux(Stmp,i-1,NX,-del) ;
		remap_flux(var,i,NX,-del) ;
#endif
	}
	}
	if(var==S) {
		if(q==0.) for(j=0;j<NY;j++) var[-3][j] = var[-3+NX][j] ;
		else {
#if defined(LINEAR)
        	        for(j=0;j<NY;j++) var[i][j] = var[i+NX][j] ;
                	remap_linear(var,i,-del) ;
#elif defined(FOURIER)
	                for(j=0;j<NY;j++) var[i][j] = var[i+NX][j] ;
	 		remap_fourier(var,i,-del) ;
#elif defined(FLUX)
			remap_flux(var,-3,NX,-del) ;
#endif
		}
	}

#ifndef SWEEPYM
        if(var==vy&&q!=0.) for(j=0;j<NY;j++) {
                var[NX][j] += -q*W*Lx ;
                var[NX+1][j] += -q*W*Lx ;
                var[-1][j] += q*W*Lx ;
                var[-2][j] += q*W*Lx ;
        }
#endif
	/* y-boundary */
	for(i=-2;i<NX+2;i++) {
		var[i][-1]   = var[i][NY-1] ;
		var[i][-2]   = var[i][NY-2] ;
		var[i][NY]   = var[i][0] ;
		var[i][NY+1] = var[i][1] ;
	}

	if(var==S) {
		var[-3][-1]  = var[-3][NY-1] ;
		var[-3][-2]   = var[-3][NY-2] ;
		var[-3][NY]   = var[-3][0] ;
		var[-3][NY+1] = var[-3][1] ;
	}
}
