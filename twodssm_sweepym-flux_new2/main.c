
#include "decs.h"
#include "defs.h"

main(int argc,char *argv[])
{
	void init(),timestep(),stepvar() ;
	void diag(int flag) ;
	double tnext,mass,mass0 ;
	int i,j ;
	int nstep = 0 ;

	if(argc < 2) {
	   fprintf(stderr,"usage: twodssm tf Lx Ly eps damp pamp g2 phi is_restart\n") ;
	   exit(0) ;
	}
	else {
		sscanf(argv[1],"%lf",&tf) ;
		sscanf(argv[2],"%lf",&Lx) ;
		sscanf(argv[3],"%lf",&Ly) ;
		sscanf(argv[4],"%lf",&eps) ;
		sscanf(argv[5],"%lf",&damp) ;
		sscanf(argv[6],"%lf",&pamp) ;
		sscanf(argv[7],"%lf",&g2) ;
		sscanf(argv[8],"%lf",&phi) ;
		sscanf(argv[9],"%d",&is_restart) ;
	}

	/* perform initializations */
	init() ;

	tnext = t ;
	mass0 = 0. ;
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) mass0 += S[i][j] ;
	mass0 /= NX*NY ;

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

		mass = 0. ;
		for(i=0;i<NX;i++)
		for(j=0;j<NY;j++) mass += S[i][j] ;
		mass /= NX*NY ;
		//fprintf(stderr,"%10.5g %10.5g\n",t,mass-mass0) ;
	}

	/* do final diagnostics */
	diag(2) ;

	fprintf(stderr,"nstep: %d\n",nstep) ;

	return(0) ;
}
