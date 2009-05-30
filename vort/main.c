
#include "decs.h"
#include "defs.h"

main(int argc,char *argv[])
{
	void init(),timestep(),stepvar() ;
	void diag(int flag) ;
	double tnext,mass,mass0 ;
	int i,j ;
	int nstep = 0 ;

	/* perform initializations */
	init() ;

	/* do initial diagnostics */
	diag(0) ;


	while(t < tf) {

		/* find timestep */
		timestep() ;

		/* step variables forward in time */
		stepvar() ;

		/* perform diagnostics */
		if(t > tnext) {
			diag(1) ;
			tnext += DTl ;
		}
		nstep++ ;

		if(nstep%100 == 0)
			fprintf(stderr,"%10.5g %d %10.5g\n",t,nstep,dt) ;
	}

	/* do final diagnostics */
	diag(2) ;

	fprintf(stderr,"nstep: %d\n",nstep) ;

	return(0) ;
}
