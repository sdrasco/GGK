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

function g = g10a(e)

%from eq. 6 of gair et al
g = (61.e0/24.e0) + (63.e0/8.e0)*e^2.e0 + (95.e0/64.e0)*e^4.e0;

end

