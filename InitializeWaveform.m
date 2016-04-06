function waveform = InitializeWaveform(S)
% 
% waveform = InitializeWaveform(S);
%
% Given input structure with the following fields:
%
%          a: dimensionless spin of larger black hole [0, 1]
%         e0: initial eccentricity (0, 1)
%         p0: initial semilatus rectum [0 inf]
%  iota0_deg: initial inclination in degrees [0 180] 
%             defined as iota = 90 - sgn(L) theta_min
%         t0: t(lambda=0) in M             [-inf inf]
%         r0: r(lambda=0) in M             [p/(1+e)  p/(1-e)]
% theta0_deg: theta(lambda=0)  in degrees [theta_min  pi-theta_min], where
%             theta_min = 90 - iota0_deg for prograde orbits (iota < 90),
%             theta_min = iota0_deg - 90 for retrograde orbits (iota > 90).
%   phi0_deg: phi(lambda=0) in degrees    [0  360] 
% sign_rdot0: sign of initial radial velocity [-1 or 1]
% sign_Tdot0: sign of initial theta velocity [-1 or 1]
%          M: mass of the large black hole, in solar masses (mu inf]
%         mu: mass of the small black hole, in solar masses [0 M)
%      tspan: final time minus initial time, in seconds (-inf inf)
% SmallSteps: number of time steps for waveform [1 inf]
%   BigSteps: number of time steps for the trajectory of the orbital
%             elements [3 inf]
%        tol: accuracy to which we will compute the orbital paramaters of 
%             the snapshots (including the coefficients in their series 
%             expansions).
%     coords: a string telling whether to use spherical (coords =
%             'spherical') or spheroidal (coords = 'spheroidal') 
%             coordinates.
% 
% this program produces a kludged calculation of the cartesian components 
% of the metric perturbation at infinity.
%
% The first 13 fields of the output structure contain the input parameters.
% The remaining fields are:
%
%    SecPerMsun: the value used for a solar mass in seconds
%       SecPerM: M in seconds
%      SecPermu: mu in seconds
%             t: coordinate time in seconds from Curt's PEIT program
%             e: eccentricity as a function of t from PEIT
%             p: semilatus rectum as a function of t from PEIT
%      iota_deg: geometric inclination (in degrees) as a functin of t
%    SplineData: output structure from call to KLUDGEDINSPIRAL
%             x: output structure from call to KLUDGEDXOFL
%             h: output structure from call to OBESRVEWAVEFORM
%        CPUsec: computation time for creating this structure, in seconds.
%
% See also PEIT, KLUDGEDINSPIRAL, KLUDGEDXOFL OBSERVEWAVEFORM
% 
% By: Steve Drasco
% 

% strip input parameters from structure
% 
% This block also functions as a check that input structure has all the
% necessary fields
a=S.a;
e0=S.e0;
p0=S.p0;
iota0_deg=S.iota0_deg;
t0=S.t0;
r0=S.r0;
theta0_deg=S.theta0_deg;
phi0_deg=S.phi0_deg;
sign_rdot0=S.sign_rdot0;
sign_Tdot0=S.sign_Tdot0;
M=S.M;
mu=S.mu;
tspan=S.tspan;
SmallSteps=S.SmallSteps;
BigSteps=S.BigSteps;
tol=S.tol;
coords=S.coords;

% start clock to keep track of computational cost
InitialCPUTime = cputime;

% make sure that inputs are reasonable
if e0 < 0 || e0 >= 1
    error('We require 0 < e0 < 1.');
end
if p0 < 0
    error('We require a positive p0.');
end
if iota0_deg > 180 || iota0_deg < 0
    error('We require 0 < iota0_deg < 180.');
end
if M < 0
    error('We require M > 0.');
end
if mu < 0 || mu > M
    error('We require 0 < mu < M.');
end
if SmallSteps < 1
    error('We require SmallSteps > 1.');
end
if BigSteps < 3 || BigSteps > SmallSteps
    error('We require 2 < BigSteps < SmallSteps.');
end
if ~strcmp(coords,'spherical') && ~strcmp(coords,'spheroidal')
    error('coords must be set to ''spherical'' or ''spheroidal''.');
end
if abs(sign_rdot0) ~= 1
    error('sign_rdot0 must be 1 or -1.');
end
if abs(sign_Tdot0) ~= 1
    error('sign_Tdot0 must be 1 or -1.');
end
if r0 < p0/(1 + e0) || r0 > p0/(1 - e0)
    error('InitializeWaveform: r0 is out of bounds.');
end
if iota0_deg < 90
  theta_min = 90 - iota0_deg;
else
  theta_min = iota0_deg - 90;
end
if theta0_deg < theta_min || theta0_deg > (180-theta_min)
    error('InitializeWaveform: theta0_deg is out of bounds.');
end

% put the input parameters into the output structure
waveform.a = a;
waveform.e0 = e0;
waveform.p0 = p0;
waveform.iota0_deg = iota0_deg;
waveform.iota_hughes_deg0 = HughesIota(iota0_deg, a, e0, p0);
waveform.t0 = t0;
waveform.r0 = r0;
waveform.theta0_deg = theta0_deg;
waveform.phi0_deg = phi0_deg;
waveform.sign_rdot0=sign_rdot0;
waveform.sign_Tdot0=sign_Tdot0;
waveform.M = M;
waveform.mu = mu;
waveform.SmallSteps = SmallSteps;
waveform.BigSteps = BigSteps;
waveform.tol = tol;
waveform.coords = coords;
waveform.SecPerMsun = 4.9255e-6;
waveform.SecPerM = waveform.SecPerMsun * waveform.M;
waveform.SecPermu = waveform.SecPerMsun * waveform.mu;

% evolve the principle orbital elements: 
% make a trajectory e(t) p(t) iota(t), using a big time step
[t, p, e, iota_hughes_rad] = peit(p0,e0,waveform.iota_hughes_deg0*pi/180,...
    t0*waveform.SecPerM,t0*waveform.SecPerM + tspan,BigSteps,M,a*M,mu);
waveform.t = t;
waveform.p = p;
waveform.e = e;

% translate to geometric definition of inclination
iota_deg = zeros(size(iota_hughes_rad));
for i=1:length(iota_hughes_rad)
    iota_deg(i)=GeometricIota(iota_hughes_rad(i)*180/pi, a, e(i), p(i), 1e-11);
end
waveform.iota_deg = iota_deg;

% translate t into units of M
t = t / waveform.SecPerM; 

% compute structure of spline data describing the evolution of all the
% orbital elements.  
waveform.SplineData = KludgedInspiral(a, M, mu, t, e, p, iota_deg, t0, tol);

% build array for the lambda values at which we'll evaulate the waveform
lambda = spline(t, waveform.SplineData.lambda,linspace(t(1),t(end),SmallSteps));

% calculate worldline (and its first three derivatives) 
% in Boyer-Lindquist coordinates.  Positional elements are evolved in this
% step.
waveform.x = KludgedXofl(waveform.SplineData,lambda,t0,...
    r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0);

% calculate cartesian components of the worldline and derivatives
if strcmp(waveform.coords,'spherical')
    % approximating Boyer-Lindquist as sphericlal flat coordinates.
    Z = CartesianWorldline_Spherical(waveform.x);
else
    % approximating Boyer-Lindquist as spheroidal flat coordinates.
    Z = CartesianWorldline_Spheroidal(waveform.x, waveform.a);
end

% put cartesian worldline into waveform structure
waveform.x.x = Z.x;
waveform.x.y = Z.y;
waveform.x.z = Z.z;
waveform.x.dx = Z.dx;
waveform.x.dy = Z.dy;
waveform.x.dz = Z.dz;
waveform.x.ddx = Z.ddx;
waveform.x.ddy = Z.ddy;
waveform.x.ddz = Z.ddz;
waveform.x.dddx = Z.dddx;
waveform.x.dddy = Z.dddy;
waveform.x.dddz = Z.dddz;
clear Z;

% calculate components of the cartesian metric at infinity
waveform.h = CartesianMetric(waveform.x);

% shift the coarse record of t so that it has the right initial value
% tshift = waveform.x.t(1) * waveform.SecPerM;
% waveform.t = waveform.t + tshift;

% log the computational cost of this job
waveform.CPUsec = cputime - InitialCPUTime;

end
