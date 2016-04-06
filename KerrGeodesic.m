function orbit = KerrGeodesic(a, e, p, iota_deg, tol)
%
% orbit = KerrGeodesic(a, e, p, iota_deg, tol)
%
% Given the dimensionless parameters
%
%         a = dimensionless black hole spin parameter ( 0 <= a < 1)
%         e = eccentricity 
%         p = dimensionless semilatus-rectum
%  iota_deg = inclination, in degrees [defined as: pi/2 - sgn(L) theta_min]
%       tol = requested fractional accuracy
%
% this program produces a structure "orbit" that contains 
% a collection of data that can describes other details of 
% the orbit.
%
% Until I improve the documentation, the fields of the orbit structure will
% hopefully be self explanatory.  An example run looks like this:
% >> orbit = KerrGeodesic(0.9, 0.5, 7, 45, 1e-4)
%
% orbit = 
%
%                  a: 0.9000
%                  e: 0.5000
%                  p: 7
%           iota_deg: 45
%    hughes_iota_deg: 45.1138
%                  E: 0.9503
%                  L: 2.2185
%                  Q: 4.9609
%               rmin: 4.6667
%               rmax: 14
%             theta_: 0.7854
%           OmegaPhi: 0.0389
%         OmegaTheta: 0.0357
%             OmegaR: 0.0243
%              Gamma: 88.2602
%         UpsilonPhi: 3.4369
%       UpsilonTheta: 3.1468
%           UpsilonR: 2.1414
%                 LP: 1.8282
%                 LT: 1.9967
%                 LR: 2.9342
%               kmax: 10
%               nmax: 10
%                 rn: [11x1 double]
%             thetak: [11x1 double]
%                 tk: [10x1 double]
%                 tn: [10x1 double]
%               phik: [10x1 double]
%               phin: [10x1 double]
%             CPUsec: 3.0597
%
% See also XOFL.
%
% Steve Drasco
% 11 June 2008
%

% Some of the comments refer to these papers:
%
% [1] W. Schmidt Class. Quantum Grav. 19 2743 (2002)
% [2] S. A. Hughes, Phys. Rev. D 61 084004 (2000)
% [3] S. Drasco and S. A. Hughes, astro-ph/0308479 
%
% WARNING: the notation conventions in these references [1-3] are NOT the 
% same. For example, z_[2,3] = ( z_[1] )^2.  This code has a mix of the 
% different conventions.

%%%%%%%%%%%% I should do this in a more object oriented way %%%%%%%%%

% start tracking CPU cost of construction
InitialCPUtime = cputime;

% define limits for indces in series expansions
nmax = 15;
kmax = 15;

% define "prograde" and theta_
iota_rad = pi * iota_deg / 180;
if iota_rad < pi/2
  prograde = 1;
  theta_ = pi/2 - iota_rad;
else
  prograde = 0;
  theta_ = iota_rad - pi/2;
end

% this routine breaks if the orbit is circular
if isequal(e, 0)
  error('This program is for non-circular orbits.');
end

% define z_ as in [1] -- NOT AS IN [2,3]
z_ = cos(theta_);

% define r1 (rmin) and r2 (rmax)
rmin = p / (1.0 + e);
rmax = p / (1.0 - e);
r1 = rmin;
r2 = rmax;

% define Schmidt's functions f(r) etc. for r1
Delta1 = r1*r1 - 2.0*r1 + a*a;
f1 = r1*r1*r1*r1 + a*a*( r1*(r1+2.0) + z_*z_*Delta1 );
g1 = 2.0*a*r1;
h1 = r1*(r1 - 2.0) + z_*z_*Delta1/(1.0-z_*z_);
d1 = (r1*r1 + a*a*z_*z_)*Delta1;

% define Schmidt's functions f(r) etc. for r2
Delta2 = r2*r2 - 2.0*r2 + a*a;
f2 = r2*r2*r2*r2 + a*a*( r2*(r2+2.0) + z_*z_*Delta2 );
g2 = 2.0*a*r2;
h2 = r2*(r2 - 2.0) + z_*z_*Delta2/(1.0-z_*z_);
d2 = (r2*r2 + a*a*z_*z_)*Delta2;

% compute determinents
kappa =          d1*h2 - d2*h1;
epsilon =        d1*g2 - d2*g1;
rho =            f1*h2 - f2*h1;
eta =            f1*g2 - f2*g1;
sigma =          g1*h2 - g2*h1;

% compute energy 
E = kappa*rho + 2.0*epsilon*sigma;
E = E +  ((-1)^prograde) *2.0* sqrt(sigma*(sigma*epsilon*epsilon ...
    + rho*epsilon*kappa - eta*kappa*kappa) );
E = E / ( rho*rho + 4.0*eta*sigma );
E = sqrt(E);

% compute Lz 
%
% for non-zero a, we can make an equation which is 
% linear in Lz from Schmidt's Eq. (B.17), by scaling
% the two equations and subtracting them so as to 
% knock out the Lz^2 term to get  
%
%  Lz = (E*E*rho - kappa) / (2.0*E*sigma).
%
% This trick doesn't work in the a=0 limit, since 
% then the "linear" equation doesn't involve Lz.
% So, to allow for a=0, we use this ugly block
% to get Lz...
%
if a > 0.001
  Lz = (E*E*rho - kappa) / (2.0*E*sigma);
else
  BB = 2.0*g1*E / h1;
  CC = (d1 - f1*E*E) / h1;
  if(prograde)
    Lz = -0.5*BB + 0.5*sqrt(BB*BB - 4.0*CC);
  else
    Lz = -0.5*BB - 0.5*sqrt(BB*BB - 4.0*CC);
  end
end

% now compute Q
Q = z_*z_*(  a*a*(1.0-E*E) + Lz*Lz/(1.0-z_*z_)  );

% inline functions for Delta and varpi^2
varpi2 = @(r) r.^2 + a^2;
Delta = @(r) r.^2 + a^2 - 2*r;

% compute d[V_r(rmin)]/dr to determine stability
dVrdr = 4*E*rmin* (E*varpi2(rmin) - a*Lz);
dVrdr = dVrdr - (2*rmin - 2)*( rmin^2 + (Lz-a*E)^2 + Q );
dVrdr = dVrdr - 2*rmin*Delta(rmin);
if dVrdr > 0
  orbit.IsStable = 1;
else
  orbit.IsStable = 0;
end

% some functions needed to compute frequencies
A = a*a*(1.0 - E*E);
B = - Lz*Lz - Q - a*a*(1.0 - E*E);
zPlusSquared = ( -B + sqrt(B*B - 4.0*A*Q) ) / (2.0*A);
zMinusSquared = z_*z_;
zmzp = zMinusSquared / zPlusSquared;
beta = a*sqrt(1.0 - E*E);
betaOvera = sqrt(1.0 - E*E);
aazPlusSquared = ( -B + sqrt(B*B - 4.0*A*Q) ) / (2.0*(1.0 - E*E));
betaRootzPlusSquared = betaOvera*sqrt(aazPlusSquared);

% some integrals in [1] needed to compute frequencies
% Note, this is the same as:
% ellE = elle(PI/2.0,sqrt(k));
% ellPi = ellpi(PI/2.0,-zMinusSquared,sqrt(k));
ellK = rf(0.0,1.0-zmzp,1.0,min(tol,1e-9));
ellE = ellK - zmzp*rd(0.0,1.0-zmzp,1.0,min(tol,1e-9))/3.0;
ellPi = ellK + zMinusSquared ...
    *rj(0.0,1.0-zmzp,1.0,1.0-zMinusSquared,min(tol,1e-9))/3.0;

% inline functions for main integrals
F = @(x) ( a*a*E/(p*p) - 2.0*a*(Lz-a*E)*(1.0 + e*cos(x))/(p*p*p) ) ...
    .*(1.0 + e*cos(x)).*(1.0 + e*cos(x)) + E;
G = @(x) Lz - 2.0*(Lz - a*E)*(1.0+e*cos(x))/p;
H = @(x) 1.0 - 2.0*(1.0 + e*cos(x))/p ...
    + a*a*(1.0 + e*cos(x)).*(1.0 + e*cos(x))/(p*p);
J = @(x) ( (1.0-E*E)*(3.0+e*e)/(1.0-e*e) - 4.0/p ...
    + (a*a*(1.0-E*E) + Lz*Lz + Q)*(1.0-e*e)/(p*p) )...
    *(1.0+e*cos(x)).*(1.0+e*cos(x)) + 2.0*(1.0 - E*E - (1.0-e*e)/p) ...
    *(1.0+e*cos(x)) + (1.0-E*E)*(1.0-e*e);
Wfunc = @(x,y) p*p*F(x)./((1.0 + e*cos(x)).*(1.0 ...
    + e*cos(x)).*H(x).*sqrt(J(x)));
Xfunc = @(x,y) 1.0./sqrt(J(x));
Yfunc = @(x,y) p*p./((1.0 + e*cos(x)).*(1.0 + e*cos(x)).*sqrt(J(x)));
Zfunc = @(x,y) G(x) ./ (H(x).*sqrt(J(x)));
% note: for the previous four functions, the second argument (y) is 
% needed if integrator is ode45, but not needed if it's quadl.

% main batch of integrals from [1]
%options = odeset('RelTol',min(tol,1e-9));
%[t,W] = ode45(Wfunc,[0 pi], 0, options); W = W(end);
%[t,X] = ode45(Xfunc,[0 pi], 0, options); X = X(end);
%[t,Y] = ode45(Yfunc,[0 pi], 0, options); Y = Y(end);
%[t,Z] = ode45(Zfunc,[0 pi], 0, options); Z = Z(end);
W = quadl(Wfunc, 0, pi, min(tol,1e-9));
X = quadl(Xfunc, 0, pi, min(tol,1e-9));
Y = quadl(Yfunc, 0, pi, min(tol,1e-9));
Z = quadl(Zfunc, 0, pi, min(tol,1e-9));

% another function needed to compute frequencies
Lambda = (Y + aazPlusSquared*X)*ellK;
Lambda = Lambda -aazPlusSquared*X*ellE;

% action-angle frequencies
OmegaR = pi*p*ellK / ((1.0 - e*e)*Lambda);
OmegaTheta = pi*betaRootzPlusSquared*X / (2.0*Lambda);
OmegaPhi = ((Z - Lz*X)*ellK + Lz*X*ellPi  ) / Lambda;

% Boyer-Lindquist frequencies
gamma = (W + aazPlusSquared*E*X)*ellK;
gamma = gamma - aazPlusSquared*E*X*ellE;
gamma = gamma / Lambda;
OmegaR = OmegaR / gamma;
OmegaTheta = OmegaTheta /gamma;
OmegaPhi = OmegaPhi / gamma;

% Mino time frequencies 
UpsilonTheta = 0.5*pi*betaRootzPlusSquared / ellK;
Gamma = UpsilonTheta / OmegaTheta;
UpsilonR = Gamma*OmegaR;
UpsilonPhi = Gamma*OmegaPhi;

% inline functions for r_n calculation
w_r = @(psi) UpsilonR*(1.0 - e*e) ...
    * vector_quadl(Xfunc, 0, psi, min(tol,1e-9)) / p;
RofPsi = @(psi) p ./ ( 1.0 + e*cos(psi) );
dwdpsi = @(psi) UpsilonR*(1.0 - e*e) ./ ( p*sqrt(J(psi)) );

% use splines for some functions in integrands that are very smooth, and
% that would otherwise be somewhat expensive to evaluate.
SplinePhase = linspace(0,pi,500);
w_r_PP = spline(SplinePhase,w_r(SplinePhase));
Spline_w_r = @(psi) ppval(w_r_PP,psi);
w_theta_PP = spline(SplinePhase,w_theta(SplinePhase,ellK,betaRootzPlusSquared,UpsilonTheta,zmzp,tol));
Spline_w_theta = @(chi) ppval(w_theta_PP,chi);

% version with spline w_r(psi)
r_n_integrand = @(psi,n) dwdpsi(psi) .* RofPsi(psi) .* cos(n*Spline_w_r(psi));

% version withOUT spline w_r(psi)
%r_n_integrand = @(psi,n) dwdpsi(psi) .* RofPsi(psi) .* cos(n*w_r(psi));

% compute coefficients r_n 
rn = zeros(nmax+1,1);
for n=0:nmax
  rn(n+1) = quadl(@(x)r_n_integrand(x,n), 0, pi, min(tol,1e-9)) / pi;
  %rn(n+1) = quadgk(@(x)r_n_integrand(x,n), 0, pi, 'RelTol', min(tol,1e-9)) / pi;
end

% inline functions for theta_k calculation
Theta = @(chi) acos( sqrt(zMinusSquared)*cos(chi) );
dwdchi = @(chi) pi ./ ( 2*ellK*sqrt(1.0 - zmzp*cos(chi).*cos(chi)) );

% version with spline w_theta(chi)
theta_k_integrand = @(chi,k) dwdchi(chi) .* Theta(chi) ...
  .* cos(k*Spline_w_theta(chi));

% version withOUT spline w_theta(chi)
%theta_k_integrand = @(chi,k) dwdchi(chi) .* Theta(chi) ...
%  .* cos(k*w_theta(chi,ellK,betaRootzPlusSquared,UpsilonTheta,zmzp,tol));

% compute coefficients theta_k
thetak = zeros(kmax+1,1);
for k=0:kmax
  thetak(k+1) = quadl(@(x)theta_k_integrand(x,k), 0, pi, min(tol,1e-9)) ...
      / pi;
  %thetak(k+1) = quadgk(@(x)theta_k_integrand(x,k), 0, pi,'RelTol', min(tol,1e-9)) ...
  %    / pi;

end

% inline functions for dt_k and dt_n calculation
dt_theta = @(theta) -a^2*E*sin(theta).^2;
dt_r = @(r) E*varpi2(r).^2 ./ Delta(r) + a*Lz*(1-varpi2(r)./Delta(r));
tk_integrand = @(w,k) ...
    dt_theta(ThetaOfLam(w/UpsilonTheta,UpsilonTheta,thetak,kmax)) ...
    .* cos(k*w);
tn_integrand = @(w,n) dt_r(RofLam(w/UpsilonR,UpsilonR,rn,nmax)) ...
    .* cos(n*w);

% compute coefficients dt_k
tk = zeros(kmax,1);
for k=1:kmax
  tk(k) = 2*quadl(@(x)tk_integrand(x,k), 0, pi, min(tol,1e-9)) ...
      / (k*pi*UpsilonTheta);
end

% compute coefficients dt_n
tn = zeros(nmax,1);
for n=1:nmax
  tn(n) = 2*quadl(@(x)tn_integrand(x,n), 0, pi, min(tol,1e-9)) ...
      / (n*pi*UpsilonR);
end

% inline functions for phi_k and phi_n calculation
dphi_theta = @(theta) Lz*csc(theta).^2;
dphi_r = @(r) a*E*(varpi2(r)./Delta(r) - 1) - a^2*Lz./Delta(r);
phik_integrand = @(w,k) ...
    dphi_theta(ThetaOfLam(w/UpsilonTheta,UpsilonTheta,thetak,kmax)) ...
    .* cos(k*w);
phin_integrand = @(w,n) dphi_r(RofLam(w/UpsilonR,UpsilonR,rn,nmax)) ...
    .* cos(n*w);

% compute coefficients phi_k
phik = zeros(kmax,1);
for k=1:kmax
  phik(k) = 2*quadl(@(x)phik_integrand(x,k), 0, pi, min(tol,1e-9)) ...
      / (k*pi*UpsilonTheta);
end

% compute coefficients dt_n
phin = zeros(nmax,1);
for n=1:nmax
  phin(n) = 2*quadl(@(x)phin_integrand(x,n), 0, pi, min(tol,1e-9)) ...
      / (n*pi*UpsilonR);
end


% fill the structure
orbit.a = a;
orbit.e = e;
orbit.p = p;
orbit.iota_deg = iota_deg;
orbit.hughes_iota_deg = acos(Lz/sqrt(Lz*Lz + Q)) * 180/pi;
orbit.E = E;
orbit.L = Lz;
orbit.Q = Q;
orbit.rmin = rmin;
orbit.rmax = rmax;
orbit.theta_ = theta_;
orbit.OmegaPhi = OmegaPhi;
orbit.OmegaTheta = OmegaTheta;
orbit.OmegaR = OmegaR;
orbit.Gamma = Gamma;
orbit.UpsilonPhi = UpsilonPhi;
orbit.UpsilonTheta = UpsilonTheta;
orbit.UpsilonR = UpsilonR;
orbit.LP = 2*pi/UpsilonPhi;
orbit.LT = 2*pi/UpsilonTheta;
orbit.LR = 2*pi/UpsilonR;
orbit.kmax = kmax;
orbit.nmax = nmax;
orbit.rn = rn;
orbit.thetak = thetak;
orbit.tk = tk;
orbit.tn = tn;
orbit.phik = phik;
orbit.phin = phin;
orbit.CPUsec = cputime - InitialCPUtime;
orbit.tol = tol;


% a version of quadl which takes a vector input 
% for the integral's outer boundary
%
% note: have turned off warnings that appear to be harmless.  
%       Should investigate/fix them instead.  
function out = vector_quadl(func, initial, vector_final, tol)

out = zeros(size(vector_final));
for i=1:length(vector_final);
  warning off;
  out(i) = quadl(func, initial, vector_final(i), tol);
  warning on;
end

% the routine for w_theta(chi) is too complicated to be inline
function out = w_theta(chi,ellK,betaRootzPlusSquared,UpsilonTheta,zmzp,tol)

% inline function called differently depending on 
% which region of the orbit we're in
lambda0 = @(x) ( ellK - ellf(0.5*pi - x, sqrt(zmzp),tol) ) ...
    / betaRootzPlusSquared;

out = zeros(size(chi));
for i=1:length(chi)

  % matlab's quadl routine will call this with arguments that are 
  % slightly out of bounds
  if(chi(i)<0)
    chi(i) = 0;
  elseif(chi(i) > 2*pi)
    chi(i) = 2*pi;
  end

  if(chi(i) >= 0 && chi(i) <= 0.5*pi)
    lambda = lambda0(chi(i));
  elseif(chi(i) > 0.5*pi && chi(i) <= pi)
    lambda = 2.0*ellK/betaRootzPlusSquared - lambda0(pi - chi(i));
  elseif(chi(i) > pi && chi(i) <= 1.5*pi)
    lambda = 2.0*ellK/betaRootzPlusSquared + lambda0(chi(i) - pi);
  elseif(chi(i) > 1.5*pi && chi(i) <= 2.0*pi)
    lambda = 4.0*ellK/betaRootzPlusSquared - lambda0(2.0*pi - chi(i));
  else
    error('chi out of range in w_theta(chi).  we require 0 < chi < 2pi.');
  end

  out(i) = UpsilonTheta*lambda;

end

% ellf eliptic integral 
function out = ellf(phi, ak, tol)

s=sin(phi);
out = s*rf(cos(phi)*cos(phi), (1.0-s*ak)*(1.0+s*ak), 1.0, min(tol,1e-9));

% theta(lambda)
function out = ThetaOfLam(lambda, UpsilonTheta, thetak, kmax)

k = (1:kmax+1)';
out = zeros(size(lambda));
for i=1:length(lambda)
  out(i) = thetak(1) ...
    + 2*sum(  thetak(2:end) .* cos((k(2:end)-1)*UpsilonTheta*lambda(i))  );
end

% r(lambda)
function out = RofLam(lambda, UpsilonR, rn, nmax) 

n = (1:nmax+1)';
out = zeros(size(lambda));
for i=1:length(lambda)
  out(i) = rn(1) ...
    + 2*sum(  rn(2:end) .* cos((n(2:end)-1)* UpsilonR*lambda(i))  );
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% eliptic integrals %%%%%%%%%%%%%
%%%%%%% not by SD, but found online %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% rj(x, y, z, p, errtol)
%
% Inputs:
%
%   x       Input vector size 1xN.
%   y       Input vector size 1xN.
%   z       Input vector size 1xN.
%   p       Input vector size 1xN.
%   errtol  Error tolerance.
%
% Matlab function to compute Carlson's symmetric elliptic integral Rj.
% Implementation of Carlson's Duplication Algorithm 3 in "Computing
% Elliptic Integrals by Duplication," by B. C. Carlson, Numer. Math.
% 33, 1-16 (1979).
%
% Returns NaN's for any argument values outside input range.
%
% Algorithm is also from Carlson's ACM TOMS Algorithm 577.
%
% This code is a complete rewrite of the algorithm in vectorized form.
% It was not produced by running a FORTRAN to Matlab converter.
%
% The following text is copied from ACM TOMS Algorithm 577 FORTRAN code:
%
%   X AND Y ARE THE VARIABLES IN THE INTEGRAL RC(X,Y).
%
%   ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
%   RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
%   16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).
%
%   SAMPLE CHOICES:  ERRTOL     RELATIVE TRUNCATION
%                               ERROR LESS THAN
%                    1.D-3      4.D-18
%                    3.D-3      3.D-15
%                    1.D-2      4.D-12
%                    3.D-2      3.D-9
%                    1.D-1      4.D-6
%
% Note by TRH:
%
%   Absolute truncation error when the integrals are order 1 quantities
%   is closer to errtol, so be careful if you want high absolute precision.
%
% Thomas R. Hoffend Jr., Ph.D.
% 3M Company
% 3M Center Bldg. 236-GC-26
% St. Paul, MN 55144
% trhoffendjr@mmm.com
%

function f = rj(x, y, z, p, errtol)

% Argument limits as set by Carlson:
LoLim = 5.0 * realmin;
UpLim = 5.0 * realmax;

% Check input arguments for acceptability:
mask = (min([x; y; z]) >= 0) & ...
       (min([(x + y); (x + z); (y + z); p]) >= LoLim) & ...
       (max([x; y; z; p]) < UpLim);

% Define internally acceptable variable ranges for iterations:
Xi = x(mask);
Yi = y(mask);
Zi = z(mask);
Pi = p(mask);

% Carlson's duplication algorithm for Rj:
Xn = Xi;
Yn = Yi;
Zn = Zi;
Pn = Pi;
sigma = 0.0;
power4 = 1.0;
etolrc = 0.5 * errtol;
Mu = (Xn + Yn + Zn + Pn + Pn) * 0.2;
Xndev = (Mu - Xn) ./ Mu;
Yndev = (Mu - Yn) ./ Mu;
Zndev = (Mu - Zn) ./ Mu;
Pndev = (Mu - Pn) ./ Mu;
epslon = max( abs([Xndev Yndev Zndev Pndev]) );
while (epslon >= errtol)
    Xnroot = sqrt(Xn);
    Ynroot = sqrt(Yn);
    Znroot = sqrt(Zn);
    lambda = Xnroot .* (Ynroot + Znroot) + Ynroot .* Znroot;
    alpha = Pn .* (Xnroot + Ynroot + Znroot) + Xnroot .* Ynroot .* Znroot;
    alpha = alpha .* alpha;
    beta = Pn .* (Pn + lambda) .* (Pn + lambda);
    sigma = sigma + power4 .* rc(alpha, beta, etolrc);
    % Here we might need to shrink the size of the arrays if alpha, beta
    % out of range for Rc.
    mask = ~isnan(sigma);
    power4 = 0.25 * power4;
    Xn = 0.25 * (Xn(mask) + lambda(mask));
    Yn = 0.25 * (Yn(mask) + lambda(mask));
    Zn = 0.25 * (Zn(mask) + lambda(mask));
    Pn = 0.25 * (Pn(mask) + lambda(mask));
    Mu = (Xn + Yn + Zn + Pn + Pn) * 0.2;
    Xndev = (Mu - Xn) ./ Mu;
    Yndev = (Mu - Yn) ./ Mu;
    Zndev = (Mu - Zn) ./ Mu;
    Pndev = (Mu - Pn) ./ Mu;
    epslon = max( abs([Xndev Yndev Zndev Pndev]) );
end
C1 = 3.0 / 14.0;
C2 = 1.0 / 3.0;
C3 = 3.0 / 22.0;
C4 = 3.0 / 26.0;
EA = Xndev .* Yndev - Zndev .* Zndev;
EB = Xndev .* Yndev .* Zndev;
EC = Pndev .* Pndev;
E2 = EA - 3.0 * EC;
E3 = EB + 2.0 * Pndev .* (EA - EC);
S1 = 1.0 + E2 .* (-C1 + 0.75D0 * C3 * E2 - 1.5D0 * C4 * E3);
S2 = EB .* (0.5D0 * C2 + Pndev .* (-C3 -C3 + Pndev*C4));
S3 = Pndev .* EA .* (C2 - Pndev*C3) - C2*Pndev .* EC;
f(mask) = 3.0 * sigma + power4 * (S1 + S2 + S3) ./ (Mu .* sqrt(Mu));

% Return NaN's where input argument was out of range:
f(~mask) = NaN;

%
% rd(x, y, z, errtol)
%
% Inputs:
%
%   x       Input vector size 1xN.
%   y       Input vector size 1xN.
%   z       Input vector size 1xN.
%   errtol  Error tolerance.
%
% Matlab function to compute Carlson's symmetric elliptic integral Rd.
% Implementation of Carlson's Duplication Algorithm 4 in "Computing
% Elliptic Integrals by Duplication," by B. C. Carlson, Numer. Math.
% 33, 1-16 (1979).
%
% Returns NaN's for any argument values outside input range.
%
% Algorithm is also from Carlson's ACM TOMS Algorithm 577.
%
% This code is a complete rewrite of the algorithm in vectorized form.
% It was not produced by running a FORTRAN to Matlab converter.
%
% The following text is copied from ACM TOMS Algorithm 577 FORTRAN code:
%
%   X AND Y ARE THE VARIABLES IN THE INTEGRAL RC(X,Y).
%
%   ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
%   RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
%   16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).
%
%   SAMPLE CHOICES:  ERRTOL     RELATIVE TRUNCATION
%                               ERROR LESS THAN
%                    1.D-3      3.D-19
%                    3.D-3      2.D-16
%                    1.D-2      3.D-13
%                    3.D-2      2.D-10
%                    1.D-1      3.D-7
%
% Note by TRH:
%
%   Absolute truncation error when the integrals are order 1 quantities
%   is closer to errtol, so be careful if you want high absolute precision.
%
% Thomas R. Hoffend Jr., Ph.D.
% 3M Company
% 3M Center Bldg. 236-GC-26
% St. Paul, MN 55144
% trhoffendjr@mmm.com
%

function f = rd(x, y, z, errtol)

% Argument limits as set by Carlson:
LoLim = 5.0 * realmin;
UpLim = 5.0 * realmax;

% Check input arguments for acceptability:
mask = (min([x; y]) >= 0) & ...
       (min([(x + y); z]) >= LoLim) & ...
       (max([x; y; z]) < UpLim);

% Define internally acceptable variable ranges for iterations:
Xi = x(mask);
Yi = y(mask);
Zi = z(mask);

% Carlson's duplication algorithm for Rf:
Xn = Xi;
Yn = Yi;
Zn = Zi;
sigma = 0.0;
power4 = 1.0;

Mu = (Xn + Yn + 3.0 * Zn) * 0.2;
Xndev = (Mu - Xn) ./ Mu;
Yndev = (Mu - Yn) ./ Mu;
Zndev = (Mu - Zn) ./ Mu;
epslon = max( abs([Xndev Yndev Zndev]) );
while (epslon >= errtol)
    Xnroot = sqrt(Xn);
    Ynroot = sqrt(Yn);
    Znroot = sqrt(Zn);
    lambda = Xnroot .* (Ynroot + Znroot) + Ynroot .* Znroot;
    sigma = sigma + power4 ./ (Znroot .* (Zn + lambda));
    power4 = 0.25 * power4;
    Xn = 0.25 * (Xn + lambda);
    Yn = 0.25 * (Yn + lambda);
    Zn = 0.25 * (Zn + lambda);
    Mu = (Xn + Yn + 3.0 * Zn) * 0.2;
    Xndev = (Mu - Xn) ./ Mu;
    Yndev = (Mu - Yn) ./ Mu;
    Zndev = (Mu - Zn) ./ Mu;
    epslon = max( abs([Xndev Yndev Zndev]) );
end
C1 = 3.0 / 14.0;
C2 = 1.0 / 6.0;
C3 = 9.0 / 22.0;
C4 = 3.0 / 26.0;
EA = Xndev .* Yndev;
EB = Zndev .* Zndev;
EC = EA - EB;
ED = EA - 6.0 * EB;
EF = ED + EC + EC;
S1 = ED .* (-C1 + 0.25 * C3 * ED - 1.50 * C4 * Zndev .* EF);
S2 = Zndev .* (C2 * EF + Zndev .* (-C3 * EC + Zndev .* C4 .* EA));
f(mask) = 3.0 * sigma + power4 * (1.0 + S1 + S2) ./ (Mu .* sqrt(Mu));

% Return NaN's where input argument was out of range:
f(~mask) = NaN;



%
% rc(x, y, errtol)
%
% Inputs:
%
%   x       Input vector size 1xN.
%   y       Input vector size 1xN.
%   errtol  Error tolerance.
%
% Matlab function to compute Carlson's symmetric elliptic integral Rc.
% Implementation of Carlson's Duplication Algorithm 2 in "Computing
% Elliptic Integrals by Duplication," by B. C. Carlson, Numer. Math.
% 33, 1-16 (1979).
%
% Returns NaN's for any argument values outside input range.
%
% Algorithm is also from Carlson's ACM TOMS Algorithm 577.
%
% This code is a complete rewrite of the algorithm in vectorized form.
% It was not produced by running a FORTRAN to Matlab converter.
%
% The following text is copied from ACM TOMS Algorithm 577 FORTRAN code:
%
%   X AND Y ARE THE VARIABLES IN THE INTEGRAL RC(X,Y).
%
%   ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
%   RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
%   16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).
%
%   SAMPLE CHOICES:  ERRTOL     RELATIVE TRUNCATION
%                               ERROR LESS THAN
%                    1.D-3      2.D-17
%                    3.D-3      2.D-14
%                    1.D-2      2.D-11
%                    3.D-2      2.D-8
%                    1.D-1      2.D-5
%
% Note by TRH:
%
%   Absolute truncation error when the integrals are order 1 quantities
%   is closer to errtol, so be careful if you want high absolute precision.
%
% Thomas R. Hoffend Jr., Ph.D.
% 3M Company
% 3M Center Bldg. 236-GC-26
% St. Paul, MN 55144
% trhoffendjr@mmm.com
%

function f = rc(x, y, errtol)

% Argument limits as set by Carlson:
LoLim = 5.0 * realmin;
UpLim = 5.0 * realmax;

% Check input arguments for acceptability:
mask = ( ((x + y) >= LoLim) & ( max([x; y]) < UpLim ) );

% Define internally acceptable variable ranges for iterations:
Xi = x(mask);
Yi = y(mask);

% Carlson's duplication algorithm:
Xn = Xi;
Yn = Yi;
Mu = (Xn + Yn + Yn) / 3.0d+0;
Sn = (Yn + Mu) ./ Mu - 2.0;
while (abs(Sn) >= errtol)
    lambda = 2.0 * sqrt(Xn) .* sqrt(Yn) + Yn;
    Xn = 0.25 * (Xn + lambda);
    Yn = 0.25 * (Yn + lambda);
    Mu = (Xn + Yn + Yn) / 3.0d+0;
    Sn = (Yn + Mu) ./ Mu - 2.0;
end
C1 = 1.0 / 7.0;
C2 = 9.0 / 22.0;
S = Sn.*Sn .* ( 0.3 + Sn .* (C1 + Sn .* (0.375D0 + Sn * C2)));
f(mask) = (1.0 + Sn) ./ sqrt(Mu);

% Return NaN's where input argument was out of range:
f(~mask) = NaN;


%
% rf(x, y, z, errtol)
%
% Inputs:
%
%   x       Input vector size 1xN.
%   y       Input vector size 1xN.
%   z       Input vector size 1xN.
%   errtol  Error tolerance.
%
% Matlab function to compute Carlson's symmetric elliptic integral Rf.
% Implementation of Carlson's Duplication Algorithm 1 in "Computing
% Elliptic Integrals by Duplication," by B. C. Carlson, Numer. Math.
% 33, 1-16 (1979).
%
% Returns NaN's for any argument values outside input range.
%
% Algorithm is also from Carlson's ACM TOMS Algorithm 577.
%
% This code is a complete rewrite of the algorithm in vectorized form.
% It was not produced by running a FORTRAN to Matlab converter.
%
% The following text is copied from ACM TOMS Algorithm 577 FORTRAN code:
%
%   X AND Y ARE THE VARIABLES IN THE INTEGRAL RC(X,Y).
%
%   ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
%   RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
%   16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).
%
%   SAMPLE CHOICES:  ERRTOL     RELATIVE TRUNCATION
%                               ERROR LESS THAN
%                    1.D-3      3.D-19
%                    3.D-3      2.D-16
%                    1.D-2      3.D-13
%                    3.D-2      2.D-10
%                    1.D-1      3.D-7
%
% Note by TRH:
%
%   Absolute truncation error when the integrals are order 1 quantities
%   is closer to errtol, so be careful if you want high absolute precision.
%
% Thomas R. Hoffend Jr., Ph.D.
% 3M Company
% 3M Center Bldg. 236-GC-26
% St. Paul, MN 55144
% trhoffendjr@mmm.com
%

function f = rf(x, y, z, errtol)

% Argument limits as set by Carlson:
LoLim = 5.0 * realmin;
UpLim = 5.0 * realmax;

% Check input arguments for acceptability:
mask = (min([x; y; z]) >= 0) & ...
       (min([(x + y); (x + z); (y + z)]) >= LoLim) & ...
       (max([x; y; z]) < UpLim);

% Define internally acceptable variable ranges for iterations:
Xi = x(mask);
Yi = y(mask);
Zi = z(mask);

% Carlson's duplication algorithm for Rf:
Xn = Xi;
Yn = Yi;
Zn = Zi;
Mu = (Xn + Yn + Zn) / 3.0d+0;
Xndev = 2.0 - (Mu + Xn) ./ Mu;
Yndev = 2.0 - (Mu + Yn) ./ Mu;
Zndev = 2.0 - (Mu + Zn) ./ Mu;
epslon = max( abs([Xndev Yndev Zndev]) );
while (epslon >= errtol)
    Xnroot = sqrt(Xn);
    Ynroot = sqrt(Yn);
    Znroot = sqrt(Zn);
    lambda = Xnroot .* (Ynroot + Znroot) + Ynroot .* Znroot;
    Xn = 0.25 * (Xn + lambda);
    Yn = 0.25 * (Yn + lambda);
    Zn = 0.25 * (Zn + lambda);
    Mu = (Xn + Yn + Zn) / 3.0d+0;
    Xndev = 2.0 - (Mu + Xn) ./ Mu;
    Yndev = 2.0 - (Mu + Yn) ./ Mu;
    Zndev = 2.0 - (Mu + Zn) ./ Mu;
    epslon = max( abs([Xndev Yndev Zndev]) );
end
C1 = 1.0 / 24.0;
C2 = 3.0 / 44.0;
C3 = 1.0 / 14.0;
E2 = Xndev .* Yndev - Zndev .* Zndev;
E3 = Xndev .* Yndev .* Zndev;
S = 1.0 + (C1 * E2 - 0.1D0 - C2 * E3) .* E2 + C3 * E3;
f(mask) = S ./ sqrt(Mu);

% Return NaN's where input argument was out of range:
f(~mask) = NaN;
