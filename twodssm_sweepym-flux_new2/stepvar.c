
#include "decs.h"

void stepvar()
{
	void step_pg() ;
	void step_visc() ;
	void step_ie() ;
	void step_trans() ;
	void damp_perturb() ;
	static int nstep = 0 ;

	/* establish numerical equilibrium before perturbing */
	damp_perturb() ;
	
	/* source steps */
	step_pg() ; 
	step_visc() ;
	if(fabs(g2 - 1.) > 1.e-6) {
		step_ie() ;
	}

	/* transport steps */
	step_trans() ; 
	
	t += dt ;
	nstep++ ;
}
