function result = Edot_mod(p,iota,e)
% Eq.20 from Gair&Glampedakis, gr-qc/0510129
%rem Matlab is case sensitive, so q and Q are different
global M spin m;
a = spin;
q = a/M;
q2 = q*q;
s = M/p;
cosi = cos(iota);
sini = sin(iota);
%for N1,N4,N5, we use Drasco's code ELzQ.m, but first using
%GeometricIota to switch from Hughes/Gair def of iota to Steve's convention
iota_deg_drasco = GeometricIota(iota*180/pi,q,e,p/M,1e-11);
[Ed,Ld,Qd] = ELzQ(q, 1.e-6, p/M, iota_deg_drasco);
E = Ed*m;
L = Ld*M*m;
Q = Qd*(m*M)^2;

N1 = E*p^4 + a*a*E*p^2 - 2*a*M*(L - a*E)*p;
N4 = (2*M*p - p^2)*L - 2*M*a*E*p;
N5 = 0.5*(2*M*p - p^2 -a^2);

% use GG's prefered equation if not very circular, otherwise go to fit.
if e > 1e-11
    result = Edot_2pn(p,iota,e) - ((1-e*e)^1.5)*(Edot_2pn(p,iota,0) +(N4/N1)*Ldot_mod(p,iota,0) ...
        +(N5/N1)*Qdot_mod(p,iota,0,Q) );
else
    result = (N4/N1)*Ldot_mod(p,iota,0)+(N5/N1)*Qdot_mod(p,iota,0,Q);
end
