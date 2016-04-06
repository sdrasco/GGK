function result = Qdot_over_sqrtQ_2pn(p,iota,e)
% Eq.56- from Gair&Glampedakis, gr-qc/0510129
%
global M spin m;
a = spin;
q = a/M;
q2 = q*q;
s = M/p;
cosi = cos(iota);
sini = sin(iota);

result = -12.8*(m*m/M)*(s^3.5)*((1-e*e)^1.5)*sini*( g9(e) - q*s^1.5*g10b(e)*cosi - s*g11(e) ...
    + pi*g12(e)*s^1.5 - g13(e)*s^2.0 + q2*(s^2.0)*(g14(e) - 5.625e0*sini*sini) );

%disp('stopping here');
%disp('stopped');
