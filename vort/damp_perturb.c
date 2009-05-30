
#include "decs.h"

void damp_perturb()
{
	int i,j ;
	double Y,Yp,X,Xp,h,kyp ;
	static int firstc = 0 ;

#if defined(LINTEST1) || defined(LINTEST2) || defined(LINTEST3) || defined(LINTEST4)
        kx0 = 2.*M_PI*m0/Lx ;
	ky0 = 2.*M_PI*n0/Ly ;
#endif

	if(t<=50./(damp*DMAX(W,1.))&&damp!=0.) {
		for(i=0;i<NX;i++)
		for(j=0;j<NY;j++) {
			vx[i][j] += -damp*DMAX(W,1.)*dt*vx[i][j] ;
		}
	}
	else {
		if(firstc==0&&pamp!=0.) {
		        for(i=0;i<NX;i++) {
				X = (i+0.5)*dx - 0.5*Lx ;
		        for(j=0;j<NY;j++) {
				vyneq[i][j] = vy[i][j] ;
				Sneq[i][j] = S[i][j] ;
				eneq[i][j] = e[i][j] ;
#if defined(PERTURB_RANC)
#if defined(BARO)
				h = h0*(1. + eps*sin(kxeq*X + phi)) ;
				cs = sqrt(g2*(gam - 1.)*h/gam) ;
		                vx[i][j] += pamp*cs*(ranc(0) - 0.5) ;
	        	        vy[i][j] *= (1. + pamp*(ranc(0) - 0.5)) ;
#elif defined(BARO2)
				h = h0*(1. + eps*sin(kxeq*X + phi)) ;
                                cs = sqrt(g2*h) ;
                                vx[i][j] += pamp*cs*(ranc(0) - 0.5) ;
                                vy[i][j] += pamp*(ranc(0) - 0.5) ;
#else
		                vx[i][j] += pamp*cs0*(ranc(0) - 0.5) ;
	        	        vy[i][j] += pamp*cs0*(ranc(0) - 0.5) ;
#endif
	                	S[i][j] *= (1. + pamp*(ranc(0) - 0.5)) ;
		                e[i][j] *= (1. + pamp*(ranc(0) - 0.5)) ;
#elif defined(PERTURB_SINE)
				kyp = 4.*M_PI/Ly ;
				Y = (j+0.5)*dy - 0.5*Ly ;
				Yp = j*dy - 0.5*Ly ;
				vx[i][j] += pamp*sin(kyp*Y) ;
				vy[i][j] *= (1. + pamp*sin(kyp*Yp)) ;
				S[i][j] *= (1. + pamp*sin(kyp*Y)) ;
				e[i][j] *= (1. + pamp*sin(kyp*Y)) ;
#elif defined(LINTEST1)
				Xp = i*dx - 0.5*Lx ;
				Yp = j*dy - 0.5*Ly ;
				Y = (j+0.5)*dy - 0.5*Ly ;
				S[i][j] = S0*(1.0 + (2.0*W/(cs0*cs0*ky0))*pamp*(q*ky0*ky0/(kx0*kx0+ky0*ky0) - 1.0)*sin(kx0*X + ky0*Y)) ;
				vx[i][j] = pamp*cos(kx0*Xp + ky0*Y) ;
				vy[i][j] = -(kx0*pamp/ky0)*cos(kx0*X + ky0*Yp) ;
				e[i][j] = u0*(1.0 + (2.0*W*g2/(cs0*cs0*ky0))*pamp*(q*ky0*ky0/(kx0*kx0+ky0*ky0) - 1.0)*sin(kx0*X + ky0*Y)) ;
#elif defined(LINTEST2)
				Xp = i*dx - 0.5*Lx ;
				Yp = j*dy - 0.5*Ly ;
				Y = (j+0.5)*dy - 0.5*Ly ;
				S[i][j] = S0*(1.0 + sr*cos(kx0*X + ky0*Y) - si*sin(kx0*X + ky0*Y)) ;
				vx[i][j] = vxr*cos(kx0*Xp + ky0*Y) - vxi*sin(kx0*Xp + ky0*Y) ;
				vy[i][j] = vyr*cos(kx0*X + ky0*Yp) - vyi*sin(kx0*X + ky0*Yp) ;
				e[i][j] = u0*(1.0 + g2*(sr*cos(kx0*X + ky0*Y) - si*sin(kx0*X + ky0*Y))) ;
#elif defined(LINTEST3)
				Xp = i*dx - 0.5*Lx ;
				Yp = j*dy - 0.5*Ly ;
				Y = (j+0.5)*dy - 0.5*Ly ;

				S[i][j] = S0 ;
				vx[i][j] = (kx0*pamp/ky0)*cos(kx0*Xp + ky0*Y) ;
				vy[i][j] = pamp*cos(kx0*X + ky0*Yp) ;
				e[i][j] = u0 ;
#elif defined(LINTEST4)
				Xp = i*dx - 0.5*Lx ;
				Yp = j*dy - 0.5*Ly ;
				Y = (j+0.5)*dy - 0.5*Ly ;
				S[i][j] *= 1.0 + sr*cos(kx0*X + ky0*Y) - si*sin(kx0*X + ky0*Y) ;
				vx[i][j] = vxr*cos(kx0*Xp + ky0*Y) - vxi*sin(kx0*Xp + ky0*Y) ;
				vy[i][j] = vyr*cos(kx0*X + ky0*Yp) - vyi*sin(kx0*X + ky0*Yp) ;
				//e[i][j] *= 1.0 + ur*cos(kx0*X + ky0*Y) - ui*sin(kx0*X + ky0*Y) ;
				e[i][j] += ur*cos(kx0*X + ky0*Y) - ui*sin(kx0*X + ky0*Y) ;
#endif
			}
			}
			if(pamp!=0.) {
				bound_var(vx,t) ;
			        bound_var(vy,t) ;
				bound_var(S,t) ;
				bound_var(e,t) ;
			}
			firstc = 1 ;
		}
		else ;
	}
}

