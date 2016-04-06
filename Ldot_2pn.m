function result = Ldot_2pn(p,iota,e)
% Eq.45 from Gair&Glampedakis, gr-qc/0510129
%
global M spin m;
a = spin;
q = a/M;
q2 = q*q;
s = M/p;
cosi = cos(iota);
sini = sin(iota);

result = -6.4*(m*m/M)*(s^3.5)*((1-e*e)^1.5)*( g9(e)*cosi + q*(s^1.5)*(g10a(e) - cosi*cosi*g10b(e)) ...
    -g11(e)*s*cosi + pi*(s^1.5)*g12(e)*cosi - g13(e)*(s^2.0)*cosi + q2*(s^2.0)*cosi*(g14(e) - (45.0/8.0)*sini*sini));
%disp('stopping here');
%disp('stopped');
