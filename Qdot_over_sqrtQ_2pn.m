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

end

function g = g9(e)

%from eq. 6 of gair et al
g = 1.e0 + (7.e0/8.e0)*e^2.e0;

end

function g = g11(e)

%from eq. 6 of gair et al
g = (1247.e0/336.e0) + (425.e0/336.e0)*e^2.e0;

end

function g = g12(e)

%from eq. 6 of gair et al
g = 4.e0 + (97.e0/8.e0)*e^2.e0;

end

function g = g13(e)

%from eq. 6 of gair et al
g = (44711.e0/9072.e0) + (302893.e0/6048.e0)*e^2.e0;

end

function g = g14(e)

%from eq. 6 of gair et al
g = (33.e0/16.e0) + (95.e0/16.e0)*e^2.e0;

end

function g = g10b(e)

%from eq. 6 of gair et al
g = (61.e0/8.e0) + (91.e0/4.e0)*e^2.e0 + (461.e0/64.e0)*e^4.e0;

end

