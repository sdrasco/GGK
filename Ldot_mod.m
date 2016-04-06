function result = Ldot_mod(p,iota,e)
% Eq.20 from Gair&Glampedakis, gr-qc/0510129
%
global M spin m;
a = spin;
q = a/M;
q2 = q*q;
s = M/p;
cosi = cos(iota);
sini = sin(iota);

result = Ldot_2pn(p,iota,e) - ((1-e*e)^1.5)*(Ldot_2pn(p,iota,0) - Ldot_fit(p,iota,e));
%if (p <8.0*M)
%   disp('stopping here');
%end
%disp('stopping here');
%disp('stopped');
