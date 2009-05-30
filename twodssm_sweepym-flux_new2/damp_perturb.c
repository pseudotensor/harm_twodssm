
#include "decs.h"

void damp_perturb()
{
	int i,j ;
	static int firstc = 0 ;

	if(t<=20./W) {
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) {
			vx[i][j] += -2.0*W*dt*vx[i][j] ;
		}
	}
	else {
		if(firstc==0) {
		        for(i=0;i<NX;i++)
		        for(j=0;j<NY;j++) {
		                cs0 = sqrt(g2*h0) ;
		                //vx[i][j] += pamp*cs0*(ranc(0) - 0.5) ;
	        	        //vy[i][j] *= (1. + pamp*cs0*(ranc(0) - 0.5)) ;
	                	S[i][j] *= (1. + pamp*(ranc(0) - 0.5)) ;
		                //e[i][j] *= (1. + pamp*(ranc(0) - 0.5)) ;
			}
			firstc = 1 ;		
		}
		else ;
	}

	//bound_var(vx,t) ;
	//bound_var(vy,t) ;
	bound_var(S,t) ;
	//bound_var(e,t) ;
}

