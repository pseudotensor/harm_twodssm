
/* function for shifting a column of data forward by some amount */

#include "decs.h"

void remap_fourier(REAL (*var)[NY+4], int i, double del)
{
	REAL tmp[NY] ;
	REAL freq,tmpj,pfac,tfac = 2.*M_PI*del ;
	void drealft(double data[], int n, int isign);
	int j ;

	for(j=0;j<NY;j++) {
		tmp[j] = var[i][j] ;
	}
	drealft(tmp-1,NY,1) ;

	for(j=2;j<NY;j+=2) {
		freq = ((REAL)j/(2.0*NY)) ;
		pfac = freq*tfac ;
		tmp[j] = (tmpj=tmp[j])*cos(pfac) - tmp[j+1]*sin(pfac) ;
		tmp[j+1] = tmpj*sin(pfac) + tmp[j+1]*cos(pfac) ;
	}
	pfac = 0.5*tfac ;

	drealft(tmp-1,NY,-1) ;
	for(j=0;j<NY;j++) var[i][j] = 2.0*tmp[j]/NY ;

        /* done! */
}

