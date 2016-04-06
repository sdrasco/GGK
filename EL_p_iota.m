function [E,L_over_M] = EL_p_iota(p,iota)
%This code is buggy, so not used in current implementation.

%We need to find the values of E,L,Q that correspond to (p,iota,e =0).
%The Q part is trivial.  We useEqs. 2.8 and 2.9 of Hughes, gr-qc/9910091
global M spin m;
M=1.e6;
%spin=8.e5;
spin=0.98*M;
m=1.0;
%test
r = 1.e6*M
E = m*(1 - 0.5*M/r);
iota_test= .3
L = m*sqrt(M*r)*cos(iota_test)
Q= (L*tan(iota_test))^2.0
rtest = rp_ra(E,L,Q);
global ptarget iotatarget;
a = spin;
ptarget = p;
iotatarget = iota;
x0 = [m*(1-0.5*M/p), m*sqrt(p/M)*cos(iota)];
options = optimset('MaxFunEvals',10^8,'TolX',10^-14,'TolFun',10^-14)
x = fsolve(@myfunc,x0,options);
E=x(1)
L = M*x(2)
Q = (L*tan(iota))^2.0
rt = rp_ra(E,L,Q)
rp = rt(1)
ra = rt(2)
disp('stopped')

function [y1,y2] = myfunc(x)
global M spin m;
global ptarget iotatarget;
a = spin;
r = ptarget;
pm = 1.0;
%pm = +1 corresponds to taking positive root in Eq.2.8 of Hughes; this is
%usually the right choice, but sometimes we may want to take pm = -1, for 
%strongly bound orbits and high spin.
E0=x(1)
L0=M*x(2)
Q0 = (L0*tan(iotatarget))^2.0;
Del = r^2.0 - 2.0*M*r + a^2.0
E1=(a^2*L0^2*(r-M)+ r*Del*Del)/( a*L0*M*(r^2.0-a^2.0) + pm*Del*sqrt(r^5.0*(r-3.0*M)+ ...
    a^4.0*r*(r+M) + a^2.0*r^2.0*(L0^2.0- 2.0*M*r + 2.0*r*r)))
Q1 = (((a^2.0+r^2.0)*E0 - a*L0)^2.0)/Del - (r^2.0 + (a*E0)^2.0 - 2.0*a*E0*L0 + L0^2.0)
y2 = (E0 - m*E1)
y1 = (Q0 - m^2*Q1)/M
%x2 = e02 + + 1.e2*(imag(ra)/M)^2;