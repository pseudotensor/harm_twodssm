
#include "decs.h"

void set_arrays()
{
	S     = (REAL (*) [NY+4])(& (S_matrix[3][2])) ;
#ifdef FLUX
	Stmp  = (REAL (*) [NY+4])(& (Stmp_matrix[3][2])) ;
#endif
	e     = (REAL (*) [NY+4])(& (e_matrix[2][2])) ;
	vx    = (REAL (*) [NY+4])(& (vx_matrix[2][2])) ;
	vy    = (REAL (*) [NY+4])(& (vy_matrix[2][2])) ;
	pot   = (REAL (*) [NY+4])(& (pot_matrix[2][2])) ;
	fl    = (REAL (*) [NY+4])(& (fl_matrix[2][2])) ;
	Sneq  = (REAL (*) [NY+4])(& (Sneq_matrix[2][2])) ;
	vyneq  = (REAL (*) [NY+4])(& (vyneq_matrix[2][2])) ;
	eneq  = (REAL (*) [NY+4])(& (eneq_matrix[2][2])) ;
	N     = (REAL (*) [NY+4])(& (N_matrix[2][2])) ;
	work1 = (REAL (*) [NY+4])(& (work1_matrix[2][2])) ;
	work2 = (REAL (*) [NY+4])(& (work2_matrix[2][2])) ;
	work3 = (REAL (*) [NY+4])(& (work3_matrix[2][2])) ;
	work4 = (REAL (*) [NY+4])(& (work4_matrix[2][2])) ;
	work5 = (REAL (*) [NY+4])(& (work5_matrix[2][2])) ;
	work6 = (REAL (*) [NY+4])(& (work6_matrix[2][2])) ;
}
