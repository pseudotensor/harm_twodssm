
/* restart subroutines */

#include "decs.h"

int restart()
{
	FILE *fp ;
	int i,j ; 
	void set_arrays() ;
	void zero_arrays() ;
	int Nx,Ny ;

	fp = fopen("restart.out","r") ;
	if(fp == NULL) {
		fprintf(stderr,"No restart file.\n") ;
		return(0) ;
	}

	fscanf(fp,"%d",&dump_cnt) ;
	fscanf(fp,"%d",&im_cnt) ;
	fscanf(fp,"%lf",&DTd) ;
	fscanf(fp,"%lf",&DTi) ;
	fscanf(fp,"%lf",&DTl) ;
	fscanf(fp,"%lf",&t) ;

	/* read in data */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
                fscanf(fp,"%f %f %f %f",
                        &(S[i][j]),
                        &(e[i][j]),
                        &(vx[i][j]),
                        &(vy[i][j])) ;
		//fprintf(stderr,"%d %d %f %f %f %f %f %f %f %c %f %f %f\n",i,j,S[i][j],e[i][j],vx[i][j],vy[i][j],tau[i][j],tcnew[i][j],Tm[i][j],rgm[i][j],ProPg[i][j],mintc[i][j],flr[i][j]) ;
	}

	fclose(fp) ;

	return(1) ;
}

void restart_dump()
{
	FILE *fp,*fopen() ;
	int i,j ; 

	fp = fopen("restart.out","w") ;
	if(fp == NULL) {
		fprintf(stderr,"error opening restart file.\n") ;
		exit(1) ;
	}

	fprintf(fp,"%d\n",NX) ;
	fprintf(fp,"%d\n",NY) ;
	fprintf(fp,"%d\n",dump_cnt) ;
	fprintf(fp,"%d\n",im_cnt) ;

	fprintf(fp,"%15.10g\n",DTd) ;
	fprintf(fp,"%15.10g\n",DTi) ;
	fprintf(fp,"%15.10g\n",DTl) ;

	fprintf(fp,"%15.10g\n",W) ;
	fprintf(fp,"%15.10g\n",q) ;
	fprintf(fp,"%15.10g\n",Lx) ;
	fprintf(fp,"%15.10g\n",Ly) ;
	fprintf(fp,"%15.10g\n",alphah) ;
	fprintf(fp,"%15.10g\n",alphag) ;
	fprintf(fp,"%15.10g\n",S0) ;
	fprintf(fp,"%15.10g\n",kxeq) ;
	fprintf(fp,"%15.10g\n",kyeq) ;
	fprintf(fp,"%15.10g\n",Qeq) ;

	fprintf(fp,"%15.10g\n",cour) ;
	fprintf(fp,"%15.10g\n",dt) ;
	fprintf(fp,"%15.10g\n",t) ;
	fprintf(fp,"%15.10g\n",nu_vnr) ;
	fprintf(fp,"%15.10g\n",nu_l) ;
	fprintf(fp,"%15.10g\n",G) ;

	fprintf(fp,"%15.10g\n",cool_fac) ;
	
	/* write out data */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		fprintf(fp,"%15.10g %15.10g %15.10g %15.10g\n",
			S[i][j],
			e[i][j],
			vx[i][j],
			vy[i][j]) ;
	}

	fclose(fp) ;
}
