
#include "decs.h"

void init()
{
	int i,j ;

	void set_arrays() ;
	void zero_arrays() ;
	void scaling() ;
	void init_spec() ;
	double ranc(int seed) ;
        int restart() ;

	double X,Y ;
	double Xp,Yp ;
	double dum ;
	double w ;
	double ksqmax,kmax ;
	double h,hp,N2,Seq ;
	FILE *ifp ;

	q = 1.5 ;
	S0 = 1. ;
	G = 0. ;
        h0 = 1.0 ;
	g2 = 1. ;
	cs = 1. ;
	Lx = 4. ;
	Ly = 8. ;
	W = 1. ;


	tf = 300. ;

	fprintf(stderr,"Lx,Ly,NX,NY,g2: %lf %lf %d %d %lf\n",Lx,Ly,NX,NY,g2) ;

	dx = Lx/(double)NX ;
	dy = Ly/(double)NY ;
	cour = 0.45 ;

        nu_vnr = 2. ;
        nu_l = 0. ;

	dump_cnt = 0 ;
	im_cnt = 0 ;

	/* set up arrays,details */
        set_arrays() ;
        zero_arrays() ;

	dt = 1.e-4 ;

	t = ranc(1) ;
        t = 0. ;
        DTd = tf/10. ;          /* dumping frequency */
        DTi = tf/1000. ;        /* image frequency */
        DTl = DTd/1000. ;       /* log frequency */

        dump_cnt = 0 ;
        im_cnt = 0 ;

        for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		X = (i+0.5)*dx - 0.5*Lx ;
		/*
		vx[i][j] = 0. + 0.5*(ranc(0) - 0.5);
		vy[i][j] = 0. + 0.5*(ranc(0) - 0.5);
		*/
		vx[i][j] = 0. ;
		vy[i][j] = 0. ;

		S[i][j] = S0 ;
		e[i][j] = 1.0 ;
		pot[i][j] = 0. ;
        }

	init_spec() ;

	/* enforce boundary conditions */
	bound_var(S,t) ;
	bound_var(e,t) ;
	bound_var(vx,t) ;
	bound_var(vy,t) ;
}
