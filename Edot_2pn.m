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

end


function g = g1(e)

%from eq. 6 of gair et al
g = 1 + (73.e0/24.e0)*e^2.e0 + (37.e0/96.e0)*e^4.e0;

end


function g = g2(e)

%from eq. 6 of gair et al
g = (73.e0/12.e0) + (823.e0/24.e0)*e^2.e0 + (949.e0/32.e0)*e^4.e0 + (491.e0/192.e0)*e^6.e0;

end


function g = g3(e)

%from eq. 6 of gair et al
g = (1247.e0/336.e0) + (9181.e0/672.e0)*e^2.e0;

end

function g = g4(e)

%from eq. 6 of gair et al
g = 4.e0 + (1375.e0/48.e0)*e^2.e0;

end

function g = g5(e)

%from eq. 6 of gair et al
g = (44711.e0/9072.e0) + (172157.e0/2592.e0)*e^2.e0;

end

function g = g6(e)

%from eq. 6 of gair et al
g = (33.e0/16.e0) + (359.e0/32.e0)*e^2.e0;

end
