function [p, e, iota] = p_e_iota(E,L,Q,M1,a,M2)
%Function for converting from E,L,Q to p, e, iota (Gair_Hughes) version.
%Units of E,L,Q are solar masses, solar masses^2, and solar masses^4.
%p (really p/M1), e, and iota are dimensionless.
%M1 is big BH mass in solar masses, M2 is smaller one in solar masses, and 
%a = (S/M), also in solar masses.
global M;
global spin; 
global m;
M = M1;
spin = a; 
m = M2;
rt = rp_ra(E,L,Q);
ra = rt(1);
rp = rt(2);
p = 2*ra*rp/(M*(ra+rp));
e = (ra-rp)/(ra+rp);
iota = atan2(sqrt(Q),abs(L));


