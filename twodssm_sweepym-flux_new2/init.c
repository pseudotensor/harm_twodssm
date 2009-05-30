
#include "decs.h"

void init()
{
	int i,j ;

	void set_arrays() ;
	void zero_arrays() ;
	void scaling() ;
	double ranc(int seed) ;
        int restart() ;

	double X,Y ;
	double Xp,Yp ;
	double dum ;
	double w ;
	double ksqmax,kmax ;
	double h,hp,N2,K,Seq ;

	q = 1.5 ;
	W = 1. ;
	S0 = 1. ;
	G = 0. ;
	h0 = 1.0 ;
	Lx *= sqrt(g2*h0)/W ;
	Ly *= sqrt(g2*h0)/W ;

	fprintf(stderr,"Lx,Ly,NX,NY,g2: %lf %lf %d %d %lf\n",Lx,Ly,NX,NY,g2) ;

	dx = Lx/(double)NX ;
	dy = Ly/(double)NY ;
	cour = 0.45 ;

        nu_vnr = 2.0 ;
        nu_l = 0. ;

	dump_cnt = 0 ;
	im_cnt = 0 ;

	/* set up arrays,details */
        set_arrays() ;
        zero_arrays() ;

	dt = 1.e-4 ;

        if(is_restart==1) restart() ;
        else {
        t = 0. ;
        DTd = tf/10. ;          /* dumping frequency */
        DTi = tf/500. ;         /* image frequency */
        DTl = DTd/1000. ;       /* log frequency */

        dump_cnt = 0 ;
        im_cnt = 0 ;

        kx0 = 2.*M_PI/Lx ;

        for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		X = (i+0.5)*dx - 0.5*Lx ;
		h = h0*(1. + eps*sin(kx0*X + phi)) ;
		hp = h0*eps*kx0*cos(kx0*X + phi) ;
		N2 = -hp*hp/(g2*h) ;
		N[i][j] = -sqrt(fabs(N2)) ;

		vx[i][j] = 0. ;
		vy[i][j] = hp/(2.*W) ;
		S[i][j] = S0 ;
		e[i][j] = h*S0/(g2 - 1.) ;

		pot[i][j] = 0. ;
        }

	Xmin = (asin((-1.+sqrt(1-eps*eps))/eps) - phi)/kx0 ;
        //Xmin = (asin(1./(kx0*kx0*h0*eps)) - phi)/kx0 ;
 	while(Xmin>Lx/2.)  Xmin -= Lx ;
 	while(Xmin<-Lx/2.) Xmin += Lx ;
	if(eps==1.) N2min = -2.*h0*kx0*kx0/g2 ;
	else N2min = -h0*eps*eps*kx0*kx0*cos(kx0*Xmin + phi)*cos(kx0*Xmin + phi)/(g2*(1. + eps*sin(kx0*Xmin + phi))) ;
        }

	/* enforce boundary conditions */
	bound_var(S,t) ;
	bound_var(e,t) ;
	bound_var(vx,t) ;
	bound_var(vy,t) ;

}
