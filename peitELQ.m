function [trange, pt, et, iotat, Et, Lt, Qt] = peitELQ(p0,e0,iota0,t0,tf,Ntimes,M1,spin1,M2)
%Function for finding p,e,iota of t given initial conditions p0,e0,iota0.
%p (really, p/M1), e and iota are all dimensionless, and iota is
%Gair_Hughes version
%Units of trange, t0, and tf are sec.
%M1 is big BH mass in solar masses, M2 is smaller one in solar masses, and 
%spin1= (S1/M1), also in solar masses.
%Ntimes is the number of timesteps, starting with t0 as 1 and ending with
%tf as N
global M;
global spin; 
global m;
%
M = M1;
spin = spin1; 
m = M2;
a_over_M = spin1/M1;
tol = 1.e-12;
hughes_iota_deg = iota0*180.0/pi;
iota_deg0 = GeometricIota(hughes_iota_deg, a_over_M, e0, p0, tol);
%
[E0, L0, Q0] = ELzQ(a_over_M, e0, p0, iota_deg0);
%now making E0,L0,Q0 dimensionful
E0 = E0*m;
L0 = L0*M*m;
Q0 = Q0*M*M*m*m;
[trange, Et, Lt, Qt] = ELQt(E0,L0,Q0,t0,tf,Ntimes,M1,spin1,M2);
%[pt, et, iotat] = p_e_iota(Et,Lt,Qt,M1,spin1,M2);
% functions should be re-written with .* (etc) instead of * to make above
% operation possible

Npts = length(trange);
pt = zeros(Npts,1);
et = zeros(Npts,1);
iotat = zeros(Npts,1);
for i = 1:Npts
   [pt(i), et(i), iotat(i)] = p_e_iota(Et(i),Lt(i),Qt(i),M1,spin1,M2); 
end  

% TEST
% to see if stability code has any fundamental flaws, we just return
% a simply function which should be reasonably easy to do numerical
% derivatives
% %SimpleFunction = @(mu,t) mu * t.^2 + t;
% SimpleFunction = @(M,t) cos(M * t * 1e-8);
% pt = SimpleFunction(M2,trange);
% et = SimpleFunction(M2,trange);
% iotat = SimpleFunction(M2,trange);
% Et = SimpleFunction(M2,trange);
% Lt = SimpleFunction(M2,trange);
% Qt = SimpleFunction(M2,trange);
%figure
%plot(trange,pt)
%figure
%plot(trange,et)
%figure
%plot(trange,iotat)
%disp('stopped');


