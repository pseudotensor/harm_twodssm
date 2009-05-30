
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
	static int tpx,tpxmin,tpy ;
	static double Xt,Yt,Xtp,Ytp ;
	double vxt,vxya,St,Sya,vyt,vyya,et,eya,Nt,N2t,Nya,pvt,Stav ;
	double vxta,vyta,eta,Sta ;
	double dPdx,dSdx ;
	double k,kx ;
	double B,C,D,E,F,G,WKB,WKB0 ;
	static int Nx = NX ;
	static int Ny = NY ;
	static int natp = 16 ;
	static int deln ;
	REAL alphat[NX] ;
	double h,hp,Seq ;
	FILE *dump_file ;
	FILE *im_file ;
	int i,j ;
	void dump(FILE *fp) ;
	void image(FILE *fp,int im_var) ;
	void restart_dump() ;
	char sweepym_config[10],remap_config[10],equilib_config[17],perturb_config[13] ;

	static double e_init ;
	static double m_init ;
	static double tdump ;
	static double timage ;
	static FILE *ener_file,*alpha_file,*time_file ;
	static double area ;
	static double period ;
	static double asin_arg,Xq0 ;

	if(call_code==0) {

		ener_file = fopen("ener.out","a") ;
		if(ener_file==NULL) {
			fprintf(stderr,"error opening energy output file\n") ;
			exit(1) ;
		}

		tdump = t - 1.e-12 ;
		timage = t - 1.e-12 ;
		area = Lx*Ly ;
		tpy = my_nint(Ny/2) ;
	}

	/* calculate energies */
	ek = 0. ;
	eth = 0. ;
	eg = 0. ;
	etot = 0. ;
	px = py = 0. ;
	mass = 0. ;
	vxt = vyt = St = et = Nt = vxya = vyya = Sya = eya = Nya = 0. ;
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

	alphah = 0. ;
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		/* Reynolds stress */
		dvx = 0.5*(vx[i][j]+vx[i+1][j]) ;
		dvy = 0.5*(vy[i][j]+vy[i][j+1]) ;
		alphah += S[i][j]*dvx*dvy*dx*dy ;
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

	fprintf(ener_file,"%10.5g\t%10.5g\t%10.5g\t%10.5g\t",t,ek/area,eth/area,eg/area) ;

	fprintf(ener_file,"%10.5g\t%10.5g\t%10.5g\t%10.5g\t",
		ekx/area,eky/area,px/area,py/area) ;

	fprintf(ener_file,"%10.5g\n",
		alphah) ;

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

	}

	/* make image of potential vorticity at regular intervals */
	if(t >= timage || call_code == 2) {
	  	sprintf(ifnam,"images/im%04d",im_cnt) ;
	  	im_file = fopen(ifnam,"w") ;

		if(im_file==NULL) {
			fprintf(stderr,"error opening image file\n") ;
			exit(2) ;
		}

		image(im_file,1) ;
		fclose(im_file) ;

		im_cnt++ ;
		timage += DTi ;
	}
}
