function result = Edot_2pn(p,iota,e)
% Eq.44 from Gair&Glampedakis, gr-qc/0510129
%
global M spin m;
a = spin;
q = a/M;
q2 = q*q;
s = M/p;
cosi = cos(iota);
sini = sin(iota);

result = -6.4*((m/M)^2.0)*(s^5.0)*((1-e*e)^1.5)*( g1(e) - q*s^1.5*g2(e)*cosi - s*g3(e) ...
    + pi*g4(e)*s^1.5 - g5(e)*s^2.0 + g6(e)*q2*s^2.0 -(527.e0/96.e0)*(q*s*sini)^2.0 );

%disp('stopping here');
%disp('stopped');
