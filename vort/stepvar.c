
#include "decs.h"

void stepvar()
{
	void step_pg() ;
	void step_visc() ;
	void step_trans() ;
	static int nstep = 0 ;

	/* establish numerical equilibrium before perturbing */
	/* source steps */
	step_pg() ; 
	step_visc() ;

	/* transport steps */
	step_trans() ; 
	
	t += dt ;
	nstep++ ;
}
