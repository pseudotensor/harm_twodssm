
/* 
   Do transport step.  Order of sweep direction is
   varied from step to step. 
*/

#include "decs.h"

void step_trans()
{
	void sweepx() ;
	void sweepy() ;
	void sweepym() ;

	sweepx() ;
	sweepy() ; 	/* sweepy comes second because of SSBC's */
	sweepym() ;

}
