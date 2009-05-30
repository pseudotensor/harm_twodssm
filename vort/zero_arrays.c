
#include "decs.h"

void zero_arrays()
{
	int i,j ;

	for(i=-2;i<=NX+1;i++)
	for(j=-2;j<=NY+1;j++) {
		S[i][j] = 0. ;
#ifdef FLUX
		Stmp[i][j] = 0. ;
#endif
		e[i][j] = 0. ;
		pot[i][j] = 0. ;
		vx[i][j] = 0. ;
		vy[i][j] = 0. ;
		fl[i][j] = 0. ;
		Sneq[i][j] = 0. ;
		vyneq[i][j] = 0. ;
		eneq[i][j] = 0. ;
		N[i][j] = 0. ;
		work1[i][j] = 0. ;
		work2[i][j] = 0. ;
		work3[i][j] = 0. ;
		work4[i][j] = 0. ;
		work5[i][j] = 0. ;
		work6[i][j] = 0. ;
	}

#ifdef FLUX
	/* third ghost-zone for calculating vx[-2][j] boundary */
	for(j=-2;j<=NY+1;j++) S[-3][j] = Stmp[-3][j] = 0. ;
#endif
}
