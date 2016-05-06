function [trange, pt, et, iotat, Et, Lt, Qt] = EvolveOrbit(p0,e0,iota0,t0,tf,Ntimes,M1,spin1,M2)
%Function for finding p,e,iota of t given initial conditions p0,e0,iota0.
%p (really, p/M1), e and iota are all dimensionless, and iota is
%Gair_Hughes version
%Units of trange, t0, and tf are sec.
%M1 is big BH mass in solar masses, M2 is smaller one in solar masses, and 
%spin1= (S1/M1), also in solar masses.
%Ntimes is the number of timesteps, starting with t0 as 1 and ending with
%tf as N
global M;
global spin; 
global m;
%
M = M1;
spin = spin1; 
m = M2;
a_over_M = spin1/M1;
tol = 1.e-12;
hughes_iota_deg = iota0*180.0/pi;
iota_deg0 = GeometricIota(hughes_iota_deg, a_over_M, e0, p0, tol);
%
[E0, L0, Q0] = ELzQ(a_over_M, e0, p0, iota_deg0);
%now making E0,L0,Q0 dimensionful
E0 = E0*m;
L0 = L0*M*m;
Q0 = Q0*M*M*m*m;

% integrate to get E(t) etc
[trange, Et, Lt, Qt] = ELQt(E0,L0,Q0,t0,tf,Ntimes,M1,spin1,M2);

% translate ELQ to pei.
[pt, et, iotat] = p_e_iota(Et,Lt,Qt,M1,spin1,M2);

end

function [p, e, iota] = p_e_iota(E,L,Q,M1,a,M2)
%Function for converting from E,L,Q to p, e, iota (Gair_Hughes) version.
%Units of E,L,Q are solar masses, solar masses^2, and solar masses^4.
%p (really p/M1), e, and iota are dimensionless.
%M1 is big BH mass in solar masses, M2 is smaller one in solar masses, and 
%a = (S/M), also in solar masses.
global M;
global spin; 
global m;
M = M1;
spin = a; 
m = M2;
[rp ra] = rp_ra(E,L,Q);
p = 2*ra.*rp./(M*(ra+rp));
e = (ra-rp)./(ra+rp);
iota = atan2(sqrt(Q),abs(L));

end

function [trange, Et, Lt, Qt] = ELQt(E0,L0,Q0,t0,tf,Ntimes,M1,a,M2)
%Function for finding E,L,Q of t given initial conditions E0,L0,Q0.
%Units of E,L,Q are solar masses, solar masses^2, and solar masses^4.
%Units of trange, t0, and tf are sec.
%M1 is big BH mass in solar masses, M2 is smaller one in solar masses, and 
%a = (S/M), also in solar masses.
%Ntimes is the number of timesteps, starting with t0 as 1 and ending with
%tf as N
global M;
global spin; 
global m;
%
M = M1;
spin = a; 
m = M2;
%My ODE uses solar masses as unit of time, so have to convert from sec.
AU = 499.00478370;
solar_mass = (AU^3.0)*(0.01720209895/86400)^2.0;
%above from Curt's unpublished "timing.pdf" 
t0 = t0/solar_mass;
tf = tf/solar_mass;

% set accuracy tollerance very high.  This ODE isn't expensive to solve, 
% and near merger, small errors in orbital parameters can push us over 
% separatrix into nonphysical orbits.
options = odeset('RelTol', 1e-12,'Events',@StopCondition);

% solve the ODEs
sol = ode45(@ELQdots,[t0 tf],[E0 L0 Q0], options);
display(sprintf('status: %d calls to flux ODEs with RelTol set to %0.1e',sol.stats.nfevals,options.RelTol));

% If wanted, restrict to specified number of geodesics (otherwise Ntimes
% aka BigSteps will be ignored)
if Ntimes
    trange = linspace(sol.x(1),sol.x(end),Ntimes);
    solution = deval(sol,trange); % uses ode45's native interpolation
    Et = solution(1,:);
    Lt = solution(2,:);
    Qt = solution(3,:);
else
    trange = sol.x;
    Et = sol.y(1,:);
    Lt = sol.y(2,:);
    Qt = sol.y(3,:);
end

%converting trange from solar mass to sec
trange = trange*solar_mass;

end

function[value,isterminal,direction] = StopCondition(t,y)
% Would be great to have a generic way to stop the integration.  For now
% just set to some minimum p.  For simple cases (circular, schwarszchild)
% the minimum p is easy to know, but for generic cases it isn't.  If
% a better stopping condition isn't determined, probably p_stop should be
% turned into an argument and controlled elsewhere.

global M;
E = y(1);
L = y(2); 
Q = y(3);
[rp ra] = rp_ra(E,L,Q);
p = 2*ra*rp/(ra+rp);
e = (ra-rp)/(ra+rp); % not used, but maybe needed at some point
p_stop = 6.02;

if p/M < p_stop
   value = 0; % we met the condition, will stop integration now
else
   value = 1; % any number other than 0 to say we haven't met condition
end
isterminal = 1; % stop the integration (0 keeps going, records meetin condition)
direction = 0; % direction of approach to stopping condition (we don't use)
end

function dy = ELQdots(t,y)
% ELQdots the rhs of the numerical kludge ODEs for E(t),L_z(t),Q(t)
% based on the formulae in GG = Gair Glampedakis, gr-qc/0510129
%
% These are modified 2PN flux rules (modified to better handle small
% eccentricty, large inclination, and to better fit Hughe's circular
% Teukolsky code).
%
global M spin m;
a = spin;
dy = [0,0,0];
E = y(1);
L = y(2); 
Q = y(3);
[rp ra] = rp_ra(E,L,Q);

% debugging: show me rp, ra, and p as time evolves
%plot(t,ra/M,'ko',t,rp/M,'kx',t,(2*ra*rp/(ra+rp))/M,'r+');
%drawnow

p = 2*ra*rp/(ra+rp); % note: p has dimensions
e = (ra-rp)/(ra+rp);
iota = atan2(sqrt(Q),abs(L));
q = a/M;
q2 = q*q;
s = M/p;
cosi = cos(iota);
sini = sin(iota);

% this block gives best agreement with what plots in GG
% Edot: GG's eq.20, with Edot_GHK replaced by eq 44 and "mod" for Ldot and Qdot
dy(1) = Edot_mod(p,iota,e);  

% Ldot: GG's eqs. 59,45, and 57
dy(2) = Ldot_mod(p,iota,e);  

% Qdot: GG's eqs. 60, 56, 57 and 58
dy(3) = Qdot_mod(p,iota,e,Q);

dy = dy';

end

function [rp ra] = rp_ra(E,L,Q)
%used to find perihelion and aphelion for any (E,L,Q), by finding roots of Eq. 1.
%Has been tested using Drasco code that makes the reverse transformation.
global M spin m;
a = spin;
rp = zeros(size(E));
ra = zeros(size(E));

%from eq. 6 of gair et al
c4 = E.^2 - m^2;
c3 = 2.0*M*m^2 * ones(size(c4));
c2 = 2.0*E.^2*a^2 - 2.0*a*L.*E - (L - a*E).^2 -Q - m^2*a^2;
c1 = 2.0*M*((L - a*E).^2 + Q);
%c0 = E^2*a^4 - 2*a^3*L*E + a^2*L^2 - a^2*((L - a*E)^2 + Q);
%note almost all terms in c0 cancel, except the last one.  all that
%remains is:
c0 = -Q*a^2;

for n = 1:length(E)

    p = [c4(n) c3(n) c2(n) c1(n) c0(n)];
    
    % find roots analytically
    r = roots(p);
    
    % zero any complex roots
    r(find(~isreal(r))) = 0;
    
    % note: Often the first root r(1) is ra (the larger one), and r(2)
    % is rp.  However, roots() doesn't return things sorted by size, so reorder
    % them by hand.  This is important near merger.
    r = flip(sort(r));
    
    % if we don't have two real roots, we faild.
    if r(1) == 0 || r(2) == 0
        error('rp_ra(): unable to find real values of rp and ra.');
    end
    
    ra(n) = r(1);
    rp(n) = r(2);
    
end

end

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
% works in Ldot and Qdot, but fails here in Edot.
if e > 1e-11
    result = Edot_2pn(p,iota,e) - ((1-e*e)^1.5)*(Edot_2pn(p,iota,0) +(N4/N1)*Ldot_mod(p,iota,0) ...
        +(N5/N1)*Qdot_mod(p,iota,0,Q) );
else
    result = (N4/N1)*Ldot_mod(p,iota,0)+(N5/N1)*Qdot_mod(p,iota,0,Q);
end

end

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

% use GG's prefered equation if not very circular, otherwise go to fit.
% wise?
if e > 1e-2
    result = Ldot_2pn(p,iota,e) - ((1-e*e)^1.5)*(Ldot_2pn(p,iota,0) - Ldot_fit(p,iota,e));
else
    result = Ldot_fit(p,iota,e);
end

end

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

function result = Ldot_fit(p,iota,e)
% Eq.57 from Gair&Glampedakis, gr-qc/0510129
%Note Ldot_fit doesn't actually have any dependence on the variable e; it's
%based on fitting to circular orbits in Kerr
global M spin m;
a = spin;
q = a/M;
q2 = q*q;
q3= q2*q;
q4 = q3*q;
s = M/p;
sqrt_s = s^0.5;
cosi = cos(iota);
sini = sin(iota);
cos2i = cosi*cosi;
cos3i = cos2i*cosi;
cos4i = cos3i*cosi;
cos5i = cos4i*cosi;
d1a = -10.742;  d1b = 28.5942; d1c = -9.07738; d2a = -1.42836; d2b = 10.7003;
d2c = -33.709; c1a = -28.1517; c1b = 60.9607; c1c = 40.9998; c2a = -0.348161;
c2b = 2.37258; c2c = -66.6584; c3a = -0.715392; c3b = 3.21593; c3c = 5.28888;
c4a = -7.61034; c4b = 128.878; c4c = -475.465; c5a = 12.2908; c5b = -113.125;
c5c = 306.119; c6a = 40.9259; c6b = -347.271; c6c = 886.503; c7a = -25.4831;
c7b = 224.227; c7c = -490.982; c8a = -9.00634; c8b = 91.1767; c8c = -297.002;
c9a = -0.645; c9b = -5.13592; c9c = 47.1982; f1a = -283.955; f1b = 736.209;
f2a = 483.266; f2b = -1325.19; f3a = -219.224; f3b = 634.499; f4a = -25.8203;
f4b = 82.078; f5a = 301.478; f5b = -904.161; f6a = -271.966; f6b = 827.319;

result = -6.4*(m*m/M)*(s^3.5)*( cosi + q*(s^1.5)*(61.e0/24.e0  - (61.e0/8.e0)*cos2i) ...
    -(1247.e0/336.e0)*s*cosi + 4.e0*pi*(s^1.5)*cosi - (44711.e0/9072.e0)*(s^2.0)*cosi + ...
    q2*(s^2.0)*cosi*(33.e0/16.e0 - (45.0/8.0)*sini*sini ) + (s^2.5)*  ...
    (q*(d1a + d1b*sqrt_s + d1c*s) +q3*(d2a + d2b*sqrt_s + d2c*s) + cosi*(c1a + c1b*sqrt_s +c1c*s) ...
     + q2*cosi*(c2a + c2b*sqrt_s + c2c*s) + q4*cosi*(c3a + c3b*sqrt_s + c3c*s) ...
     +q*cos2i*(c4a + c4b*sqrt_s + c4c*s) + q3*cos2i*(c5a + c5b*sqrt_s + c5c*s) ...
     +q2*cos3i*(c6a + c6b*sqrt_s + c6c*s) + q4*cos3i*(c7a + c7b*sqrt_s + c7c*s) ...
     +q3*cos4i*(c8a + c8b*sqrt_s + c8c*s) + q4*cos5i*(c9a + c9b*sqrt_s + c9c*s) ) ...
     +(s^3.5)*q*cosi*( f1a + f1b*sqrt_s + q*(f2a + f2b*sqrt_s) + q2*(f3a + f3b*sqrt_s) ...
     +cos2i*(f4a + f4b*sqrt_s) + q*cosi*cosi*(f5a + f5b*sqrt_s) + q2*cos2i*(f6a + f6b*sqrt_s)));
     
end

function result = Qdot_mod(p,iota,e,Q)
% Eq.60 from Gair&Glampedakis, gr-qc/0510129
%rem Matlab is case sensitive, so q and Q are different
global M spin m;
a = spin;
q = a/M;
q2 = q*q;
s = M/p;
cosi = cos(iota);
sini = sin(iota);

% use GG's prefered equation if not very circular, otherwise go to fit.
% wise?
if e > 1e-2
    result = sqrt(Q)*(Qdot_over_sqrtQ_2pn(p,iota,e) - ((1-e*e)^1.5)*Qdot_over_sqrtQ_2pn(p,iota,0) ...
        +2.0*((1-e*e)^1.5)*tan(iota)*(Ldot_fit(p,iota,e) + sqrtQ_oversin2i_iotadot_fit(p,iota,e) ));
else
    result = sqrt(Q)*2.0*((1-e*e)^1.5)*tan(iota)*(Ldot_fit(p,iota,e) + sqrtQ_oversin2i_iotadot_fit(p,iota,e) );
end

end

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

function result = sqrtQ_oversin2i_iotadot_fit(p,iota,e)
% Eq.58 from Gair&Glampedakis, gr-qc/0510129, multiplied by sqrt{Q}/sin^2i
%
global M spin m;
a = spin;
q = a/M;
q2 = q*q;
q3= q2*q;
q4 = q3*q;
s = M/p;
sqrt_s = s^0.5;
cosi = cos(iota);
sini = sin(iota);
cos2i = cosi*cosi;
cos3i = cos2i*cosi;
cos4i = cos3i*cosi;
cos5i = cos4i*cosi;
d1a = -10.742;  d1b = 28.5942; d1c = -9.07738; d2a = -1.42836; d2b = 10.7003;
d2c = -33.709; c1a = -28.1517; c1b = 60.9607; c1c = 40.9998; c2a = -0.348161;
c2b = 2.37258; c2c = -66.6584; c3a = -0.715392; c3b = 3.21593; c3c = 5.28888;
c4a = -7.61034; c4b = 128.878; c4c = -475.465; c5a = 12.2908; c5b = -113.125;
c5c = 306.119; c6a = 40.9259; c6b = -347.271; c6c = 886.503; c7a = -25.4831;
c7b = 224.227; c7c = -490.982; c8a = -9.00634; c8b = 91.1767; c8c = -297.002;
c9a = -0.645; c9b = -5.13592; c9c = 47.1982; f1a = -283.955; f1b = 736.209;
f2a = 483.266; f2b = -1325.19; f3a = -219.224; f3b = 634.499; f4a = -25.8203;
f4b = 82.078; f5a = 301.478; f5b = -904.161; f6a = -271.966; f6b = 827.319;

c10a = -0.0309341; c10b = -22.2416; c10c = 7.55265; c11a = - 3.33476; c11b = 22.7013;
c11c = -12.47; f7a = -162.268; f7b = 247.168; f8a = 152.125; f8b = -182.165;
f9a = 184.465; f9b = -267.553; f10a = -188.132; f10b = 254.067;


result = 6.4*(m*m/M)*q*(s^5.0)*( (61.e0/24.e0) + s*(d1a + d1b*sqrt_s + d1c*s) +q2*s*(d2a + d2b*sqrt_s + d2c*s) ...
    + q*cosi*sqrt_s*(c10a + c10b*s + c10c*s^1.5) + q2*cos2i*s*(c11a + c11b*sqrt_s + c11c*s) ...
    + q3*cosi*(s^2.5)*(f7a + f7b*sqrt_s + q*(f8a +f8b*sqrt_s) + cos2i*(f9a + f9b*sqrt_s) +q*cos2i*(f10a + f10b*sqrt_s)));

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

function g = g9(e)

%from eq. 6 of gair et al
g = 1.e0 + (7.e0/8.e0)*e^2.e0;

end

function g = g10a(e)

%from eq. 6 of gair et al
g = (61.e0/24.e0) + (63.e0/8.e0)*e^2.e0 + (95.e0/64.e0)*e^4.e0;

end

function g = g10b(e)

%from eq. 6 of gair et al
g = (61.e0/8.e0) + (91.e0/4.e0)*e^2.e0 + (461.e0/64.e0)*e^4.e0;

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

