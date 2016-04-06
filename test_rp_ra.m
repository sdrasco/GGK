function [p_curt,iota_curt,e_curt] = test_rp_ra(p,iota,e)
% This code simply tests whether Curt's codes for converting between
% (E,L,Q) and (p,iota,e) are consitent with Steve's.  We start with p,
% iota_Gair = iota_hughes, and e, and use Steve's codes to convert to
% E,L,Q.  Then we use my codes to convert back, and see if we get same
% answer.
global M spin m
M= 1.e6;
m = 1.e1;
spin = 0.5*M;
a = spin;

%for N1,N4,N5, we use Drasco's code ELzQ.m to get E,L for p,iota.
%Unfortunately, Steve's definition of iota is slightly diff from Gair's, so
%this should be improved
iota_deg_drasco = GeometricIota(iota*180.0/pi,a/M,e,p/M,10^-11);
[Ed,Ld,Qd] = ELzQ(a/M, e, p/M, iota_deg_drasco);
E = Ed*m;
L = Ld*M*m;
Q = Qd*(m*M)^2.0
rt = rp_ra(E,L,Q)
ra = rt(1);
rp = rt(2);
p_curt = 2*ra*rp/(ra+rp)
e_curt = (ra-rp)/(ra+rp)
iota_curt = atan2(sqrt(Q),abs(L))
%disp('stopping here');
%disp('stopped');
