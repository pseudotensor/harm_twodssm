h0 = 1.
gam = 3.
g2 = 2.
eps = 0.03
K = 1.
phi = Pi/2.
Lx = 2.
kx0 = 2.*Pi/Lx
q = 1.5

h = h0*(1. + eps*Sin[kx0*x + phi])
hp = h0*eps*kx0*Cos[kx0*x + phi]
hpp = -h0*eps*kx0^2.*Sin[kx0*x + phi]
s = (h*(gam - 1.)/(gam*K))^(1./(gam - 1.))
p = K*s^gam
vy = -q*x + hp/2.
N2 = -(gam - g2)*hp*hp/(g2*(gam-1.)*h)
K2 = 2.*(2. - q + hpp/2.)

Plot[N2+K2,{x,-Lx/2.,Lx/2.}]
Plot[K2,{x,-Lx/2.,Lx/2.}]
Plot[N2,{x,-Lx/2.,Lx/2.}]
Plot[{N2,K2},{x,-Lx/2.,Lx/2.}]