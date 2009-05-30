
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
		remap_linear(var,i,del) ;
	}
	for(i=-2;i<=-1;i++) {
		remap_linear(var,i,-del) ;
	}
	if(var==S) remap_linear(var,-3,-del) ;

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
