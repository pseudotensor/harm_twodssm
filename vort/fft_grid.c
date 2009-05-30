
#include "decs.h"

/* fft on grid */
void fft_grid(inr,ini,outr,outi,nx,ny,dir)
float **inr,**ini,**outr,**outi ;
int nx,ny,dir ;
{
	static float *inr1,*ini1,*outr1,*outi1 ;
	static int firstc = 1 ;
	float *vector() ;
	int i,j,n ;

	/* hold onto allocated arrays */
	if(firstc) {
		firstc = 0 ;
		if (nx > ny) n = nx ;
		else n = ny ;

		inr1 = vector(0,n-1) ;
		ini1 = vector(0,n-1) ;
		outr1 = vector(0,n-1) ;
		outi1 = vector(0,n-1) ;
	}

	/* sweep in x direction */
	for(i=0;i<NX;i++) {
		for(j=0;j<NY;j++) {
			inr1[j] = inr[i][j] ;
			ini1[j] = ini[i][j] ;
		}
		tfft(inr1,ini1,outr1,outi1,NY,dir) ;
		for(j=0;j<NY;j++) {
			outr[i][j] = outr1[j];
			outi[i][j] = outi1[j];
		}
	}

	/* sweep in y direction */
	for(j=0;j<NY;j++) {
		for(i=0;i<NX;i++) {
			inr1[i] = outr[i][j] ;
			ini1[i] = outi[i][j] ;
		}
		tfft(inr1,ini1,outr1,outi1,NX,dir) ;
		for(i=0;i<NX;i++) {
			outr[i][j] = outr1[i];
			outi[i][j] = outi1[i];
		}
	}

	/* done! */
}

