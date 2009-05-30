F[r] = S[r]*W/(K*K-w*w+N*N)
eq1 = I*w*dS[r]-(1./r)*(r*dvr[r]*S[r])'-I*kp*S[r]*dvp[r]
eq2 = -dS[r]+dh[r]*S[r]/(c*c) + I*S[r]*dvr[r]/(w*Ls[r])
eq3 = S[r]*dvr[r] - I*F[r]*(w*(dh[r]'-dh[r]/Ls[r])/W-2.*kp*dh[r])
eq4 = S[r]*dvp[r] - F[r]*(-kp*(w/W - N[r]*N[r]/(w*W))*dh[r] + K*K*(dh[r]'-dh[r]/Ls)/(2.*W*W))
Eliminate[{eq1==0.,eq2==0.,eq3==0.,eq4==0.}dS[r]]
