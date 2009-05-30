phi=-Pi/2.

eps=0.5

h=1.

L=1.1*Sqrt[2.]*Pi

K=1.

gam=3.1

Gam=3.

vy=-1.5*x+h*eps*Pi*Cos[2.*Pi*x/L+phi]/L

N2=(gam-Gam)/(gam*(Gam-1.))*(h*(2.*Pi*eps/L)^2.*Cos[2.*Pi*x/L+phi]/(1.+eps*Sin[2.*Pi*x/L+phi]))

Love=K^((gam-2.)/(gam*(Gam-1.)))*((Gam-1.)/Gam)^((2.*Gam-gam)/(gam*(Gam-1.)))*(h*(1.+eps*Sin[2.*Pi*x/L+phi]))^((2.*Gam-gam)/(gam*(Gam-1.)))/(1.-h*eps*(2.*Pi/L)^2.*Sin[2.*Pi*x/L+phi])

Lp=(1.-h*eps*(2.*Pi/L)^2.*Sin[2.*Pi*x/L+phi])*((2.*Gam-gam)/(Gam-1.))/(1.+eps*Sin[2.*Pi*x/L+phi])+h*eps*(2.*Pi/L)^2.

Plot[vy,{x,-L/2.,L/2.}]

Plot[N2,{x,-L/2.,L/2.}]

Plot[Lp,{x,-L/2.,L/2.}]

Plot[{vy,N2,Love},{x,-L/2.,L/2.}]

Plot[Love,{x,-L/2.,L/2.}]

A = m