#include "decs.h"

void dump(FILE *fp)
{
	int i,j ;
	double p,X,Save,qeff,dPdx,dSdx,N2,dvx,dvy,alpha ;
	int nghost=NGHOST ;

	fprintf(fp,"%d %d %lf %lf %lf %lf %lf\n",NX,NY,Lx,Ly,W,q,t) ;

        fprintf(fp,"#%8s %8s %10s %8s %8s %8s %8s %8s %8s\n",
		"x","y","s","p","vx","vy","pot","vort","divv") ;

        for(i=-nghost;i<NX+nghost;i++)
        for(j=-nghost;j<NY+nghost;j++) {
		if(fabs(g2 - 1.) > 1.e-6) p = (g2 - 1.)*e[i][j] ;
		else p = cs*cs*S[i][j] ;

		fprintf(fp,"%8.4g %8.4g ", (i-NX/2+0.5)*dx, 
						(j-NY/2+0.5)*dy) ;
		fprintf(fp,"%10.6g %8.4g ", S[i][j],p) ;
		fprintf(fp,"%8.4g %8.4g ", vx[i][j],vy[i][j]) ;
		fprintf(fp,"%8.4g ",pot[i][j]) ;

		/* potential vorticity */
		Save = 0.25*(S[i][j] + S[i-1][j] + S[i][j-1] + S[i-1][j-1]) ;
		fprintf(fp,"%8.4g ",((2.-q)*W + (vy[i][j] - vy[i-1][j])/dx
			- (vx[i][j] - vx[i][j-1])/dy)/Save
			) ;

		/* divv */
		fprintf(fp,"%8.4g ",(vx[i+1][j] - vx[i][j])/dx
			- (vy[i][j+1] - vy[i][j])/dy
			) ;

		/* qeff */
		qeff = q-(vy[i+1][j]-vy[i-1][j])/(2.*dx) ;
		fprintf(fp,"%8.4g ",qeff) ;

		/* N^2 */
                dPdx = (g2 - 1.)*(e[i+1][j]-e[i-1][j])/(2.*dx) ;
                dSdx = (S[i+1][j]-S[i-1][j])/(2.*dx) ;
                N2 = -(dPdx/(g2*S[i][j]))*(dPdx/((g2-1.)*e[i][j])-g2*dSdx/S[i][j]) ;
		if(i<-2||i>NX+1) N2 = 0.0 ;
                fprintf(fp,"%8.4g ",N2) ;
                
		/* alpha */
		dvx = 0.5*(vx[i][j]+vx[i+1][j]) ;
		dvy = 0.5*(vy[i][j]+vy[i][j+1]-vyneq[i][j]-vyneq[i][j+1]) ;
		alpha = S[i][j]*dvx*dvy ;
		fprintf(fp,"%8.4g ",alpha) ;
		fprintf(fp,"\n") ;
		
	}
}
