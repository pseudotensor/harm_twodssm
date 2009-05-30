
/* function for shifting a column of data forward by some amount */

#include "decs.h"

void remap_linear(REAL (*var)[NY+4], int i, double del)
{
	double tmp[NY] ;
	int j,idel,jp ;
	double fdel ;

	/* separate into integral & fraction part of shift */
        idel = (int)floor(del) ;
        fdel = del - (double)idel ;

	/* do integral part of shift */
        iremap(var[i],idel) ;

        /* do fractional part of shift */
	for(j=0;j<NY;j++) tmp[j] = var[i][j] ;
        for(j=0;j<NY;j++) {
                jp = (j+1)%NY ;
                var[i][jp] = fdel*tmp[j] + (1. - fdel)*tmp[jp] ;
        }
}


/* void iremap(REAL *var, int idel)
{
        static REAL tmp[NY] ;
	int k,kp ;

	for(k=0;k<NY;k++) tmp[k] = var[k] ;

        for(k=0;k<NY;k++) {
                kp = k + idel ;
                while(kp < 0) kp += NY ;
                while(kp > NY-1) kp -= NY ;
                var[kp] = tmp[k] ;
	}
} */
