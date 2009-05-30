
#include "decs.h"

/* diagnostics subroutine */

/*
	call codes:
	0: initial
	1: logfile call
	2: dumpfile call
*/

void diag(int call_code)
{
	char dfnam[10] ;
	char ifnam[100] ;
	double ek,eki,ekx,eky,eth,ethi,etot,X,Y,px,pxi,py,pyi,mass ;
	double dvx,dvy,em ;
	double phase ;
	double e_fin ;
	double eg,egi ;
	double m_fin ;
	double vxa,vya ;
	double Nsq ;
	static int tpy,tpx ;
	double Xt ;
	double Nx,Ny ;
	FILE *dump_file ;
	FILE *im_file ;
	int i,j ;
	void dump(FILE *fp) ;
	void image(FILE *fp,int im_var) ;
	void restart_dump() ;

	static double e_init ;
	static double m_init ;
	static double tdump ;
	static double timage ;
	static FILE *ener_file ;
	static double area ;
	static double period ;

	if(call_code==0) {
		ener_file = fopen("ener.out","a") ;
		if(ener_file==NULL) {
			fprintf(stderr,"error opening energy output file\n") ;
			exit(1) ;
		}
		tdump = t - 1.e-12 ;
		timage = t - 1.e-12 ;
		area = Lx*Ly ;
		Nx = NX ;
		Ny = NY ;
		tpx = my_nint(Nx*(Xmin/Lx + 0.5)) ;
		tpy = my_nint(Ny/2.) ;
		Xt = (tpx+0.5)*dx - 0.5*Lx ;
		if(is_restart!=1) fprintf(ener_file,"#config: Lx = %.5g, Ly = %.5g, NX = %d, NY = %d, N2min = %.5g, Xmin = %.5g, pamp = %.5g, phi = %.5g, eps = %.5g, g2 = %g, h0 = %g, G = %g\n",Lx,Ly,NX,NY,N2min,Xmin,pamp,phi,eps,g2,h0,G) ;
	}

	/* calculate energies */
	ek = 0. ;
	eth = 0. ;
	eg = 0. ;
	etot = 0. ;
	px = py = 0. ;
	mass = 0. ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		vxa = 0.5*(vx[i+1][j]+vx[i][j]) ;
		vya = 0.5*(vy[i][j+1]+vy[i][j]) ;
		pxi = vxa*S[i][j] ;
		pyi = vya*S[i][j] ;
		
		ekx = 0.5*S[i][j]*vxa*vxa ;
		eky = 0.5*S[i][j]*vya*vya ;
		eki = ekx + eky ;

		egi = 0.5*S[i][j]*pot[i][j] ;

		if(fabs(g2 - 1.) > 1.e-6) ethi = e[i][j] ;
		else ethi = cs*cs*S[i][j] ;

		ek += eki*dx*dy ;
		eth += ethi*dx*dy ;
		eg += egi*dx*dy ;
		px += pxi*dx*dy ;
		py += pyi*dx*dy ;
		mass += S[i][j]*dx*dy ;
	
	}
	etot = ek+eth+eg ;

	/* usually, alphah is calculated in potsolv.
  	   but potsolv is not called if G < 0; calculate
 	   Reynolds stress here in that case */
	if(G <= 0.) {
                alphah = 0. ;
		em = 0. ;

                for(i=NX/4;i<3*NX/4;i++)
                for(j=0;j<NY;j++) {

                	/* Reynolds stress */
                	X = (i+0.5)*dx - 0.5*Lx ;
                	dvx = 0.5*(vx[i][j]+vx[i+1][j]) ;
                	dvy = 0.5*(vy[i][j]+vy[i][j+1]) ;
                        alphah += S[i][j]*dvx*dvy*dx*dy ;
			em += e[i][j]*dx*dy ;
                }
                alphah /= g2*(g2-1.)*em*q ;
	}

	/* calculate energies */
	if(call_code == 0) {
		e_init = etot ;
		m_init = 0. ;
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) 
			m_init += S[i][j]*dx*dy ;
	}
	else if (call_code == 2) {
		e_fin = etot ;
		m_fin = 0. ;
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) 
			m_fin += S[i][j]*dx*dy ;
		fprintf(stdout,"\n\nEnergy: ini,fin,del: %g %g %g\n",
			e_init,e_fin,(e_fin-e_init)/(fabs(e_init) +
				SMALL)) ;
		fprintf(stdout,"mass: ini,fin,del: %g %g %g\n",
			m_init,m_fin,(m_fin-m_init)/(fabs(m_init) +
				SMALL)) ;
	}
	else {
#if 0
		fprintf(ener_file,"%15.15g\t%15.15g\t%15.15g\t%15.15g\t%15.15g\n",t,
			(0.25*( S[NX/2+tpx][NY/2+tpy]
			+ S[NX/2-1+tpx][NY/2+tpy]
			+ S[NX/2+tpx][NY/2-1+tpy]
			+ S[NX/2-1+tpx][NY/2-1+tpy])-S0)*1.,
			0.5*(vx[NX/2+tpx][NY/2+tpy] + vx[NX/2+tpx][NY/2-1+tpy])*1.,
			0.5*(vy[NX/2+tpx][NY/2+tpy] + vy[NX/2-1+tpx][NY/2+tpy])*1.,
			(0.25*( e[NX/2+tpx][NY/2+tpy]
			+ e[NX/2-1+tpx][NY/2+tpy]
			+ e[NX/2+tpx][NY/2-1+tpy]
			+ e[NX/2-1+tpx][NY/2-1+tpy])-u0)*1.) ;
#endif

#if 0
                fprintf(ener_file,"%15.15g\t%15.15g\t%15.15g\t%15.15g\t%15.15g\n",t,
                        S[tpx][tpy],
                        0.5*(vx[tpx][tpy] + vx[tpx+1][tpy]),
                        0.5*(vy[tpx][tpy] + vy[tpx][tpy+1]),
                        e[tpx][tpy]) ;
#endif

#if 1 
	fprintf(ener_file,"%10.5g\t%10.5g\t%10.5g\t%10.5g\t",t,ek/area,eth/area,eg/area) ;
	fprintf(ener_file,"%10.5g\t%10.5g\t%10.5g\t%10.5g\t",
		ekx/area,eky/area,px/area,py/area) ;
	fprintf(ener_file,"%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\t%10.5g\n",
		alphag,alphah,vx[tpx][tpy],S[tpx][tpy]-S0,vy[tpx][tpy]-h0*eps*kx0*cos(kx0*Xt + phi)/(2.*W),e[tpx][tpy]-h0*(1. + eps*sin(kx0*Xt + phi))*S0/(g2 - 1.)) ;
#endif
		fflush(ener_file) ;
	}

	/* dump at regular intervals */
	if(t >= tdump || call_code == 2) {
		fprintf(stderr,"t: %lf %lf %lf\n",t,tdump,DTd) ;
		sprintf(dfnam,"dump%03d",dump_cnt) ;
		dump_file = fopen(dfnam,"w") ;

		if(dump_file==NULL) {
			fprintf(stderr,"error opening dump file\n") ;
			exit(2) ;
		}

		dump(dump_file) ;
		fclose(dump_file) ;

		dump_cnt++ ;
		tdump += DTd ;

		/*
		restart_dump() ;
		*/
	}

	/* make image of density at regular intervals */
	if(t >= timage || call_code == 2) {
	  	sprintf(ifnam,"imrho%04d",im_cnt) ;
	  	im_file = fopen(ifnam,"w") ;

		if(im_file==NULL) {
			fprintf(stderr,"error opening image file\n") ;
			exit(2) ;
		}

		image(im_file,1) ;
		fclose(im_file) ;

	  	sprintf(ifnam,"imeth%04d",im_cnt) ;
	  	im_file = fopen(ifnam,"w") ;

		if(im_file==NULL) {
			fprintf(stderr,"error opening image file\n") ;
			exit(2) ;
		}

		image(im_file,2) ;
		fclose(im_file) ;

		sprintf(ifnam,"imv_x%04d",im_cnt) ;
		im_file = fopen(ifnam,"w") ;

		if(im_file==NULL) {
			fprintf(stderr,"error opening image file\n") ;
			exit(2) ;
                }
		
		image(im_file,4) ;
		fclose(im_file) ;

                sprintf(ifnam,"imv_y%04d",im_cnt) ;
                im_file = fopen(ifnam,"w") ;

                if(im_file==NULL) {
                        fprintf(stderr,"error opening image file\n") ;
                        exit(2) ;
                }

                image(im_file,5) ;
                fclose(im_file) ;

		im_cnt++ ;
		timage += DTi ;
	}
}
