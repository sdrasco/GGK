function waveform = phase_InitializeWaveform(a, e0, p0, iota_deg0, phase0, M, mu, tspan, SmallSteps, BigSteps, tol, coords)
% 
% waveform = InitializeWaveform(a, e0, p0, iota_deg0, phase0, M1, M2, tspan, SmallSteps, BigSteps, tol, coords)
%
% Given input parameters, 
%
%          a: dimensionless spin of larger black hole [0, 1]
%         e0: initial eccentricity (0, 1)
%         p0: initial semilatus rectum [0 inf]
%  iota_deg0: initial inclination in degrees [0 180] 
%             defined as iota = 90 - sgn(L) theta_min
%     phase0: 4-vector (t, r, theta, phi) of initial phases [0 2pi]
%          M: mass of the large black hole, in solar masses (mu inf]
%         mu: mass of the small black hole, in solar masses [0 M)
%      tspan: final time minus initial time, in seconds (0 inf]
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

% start clock to keep track of computational cost
InitialCPUTime = cputime;

% make sure that inputs are reasonable
if e0 < 0 || e0 >= 1
    error('We require 0 < e0 < 1.');
end
if p0 < 0
    error('We require a positive p0.');
end
if iota_deg0 > 180 || iota_deg0 < 0
    error('We require 0 < iota_deg0 < 180.');
end
%for i=1:4
%    if phase0(i) < 0 || phase0(i) > 2*pi
%        error('Though it''s not necessary, we require 0 < phase0 < 2pi.');
%    end
%end
if M < 0
    error('We require M > 0.');
end
if mu < 0 || mu > M
    error('We require 0 < mu < M.');
end
if tspan <= 0
    error('We require tspan > 0.');
end
if SmallSteps < 1
    error('We require SmallSteps > 1.');
end
if BigSteps && (BigSteps < 3 || BigSteps > SmallSteps)
    error('We require 2 < BigSteps < SmallSteps.');
end
if ~strcmp(coords,'spherical') && ~strcmp(coords,'spheroidal')
    error('coords must be set to ''spherical'' or ''spheroidal''.');
end

% put the input parameters into the output structure
waveform.a = a;
waveform.e0 = e0;
waveform.p0 = p0;
waveform.iota_deg0 = iota_deg0;
waveform.iota_hughes_deg0 = HughesIota(iota_deg0, a, e0, p0);
waveform.phase0 = phase0;
waveform.M = M;
waveform.mu = mu;
waveform.SmallSteps = SmallSteps;
waveform.BigSteps = BigSteps;
waveform.tol = tol;
waveform.coords = coords;

% make a trajectory e(t) p(t) iota(t), using a big time step
[t, p, e, iota_hughes_rad] = peit(p0,e0,waveform.iota_hughes_deg0*pi/180,0,tspan,BigSteps,M,a*M,mu);
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
waveform.SecPerMsun = 4.9255e-6;
waveform.SecPerM = waveform.SecPerMsun * waveform.M;
waveform.SecPermu = waveform.SecPerMsun * waveform.mu;
t = t / waveform.SecPerM;

% compute structure of spline data describing the evolution of all the
% orbital parameters.  
waveform.SplineData = phase_KludgedInspiral(a, M, mu, t, e, p, iota_deg, phase0, tol);

% build array for the lambda values at which we'll evaulate the waveform
% lambda = linspace(waveform.SplineData.lambda(1),waveform.SplineData.lambda(end),SmallSteps);
lambda = spline(t, waveform.SplineData.lambda, linspace(t(1),t(end),SmallSteps));

% from initial phases, compute initial lambda_{t, r, theta, phi}
MinoFrequencies0 = [waveform.SplineData.geodesic{1}.Gamma ...
    waveform.SplineData.geodesic{1}.UpsilonR ...
    waveform.SplineData.geodesic{1}.UpsilonTheta ...
    waveform.SplineData.geodesic{1}.UpsilonPhi];
lambda_x0 = phase0 ./ MinoFrequencies0;

% calculate worldline (and its first three derivatives) 
% in Boyer-Lindquist coordinates
waveform.x = phase_KludgedXofl(waveform.SplineData,lambda,lambda_x0);

% calculate cartesian compoents of the worldline and derivatives
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
t0 = waveform.x.t(1) * waveform.SecPerM;
waveform.t = waveform.t + t0;

% log the computational cost of this job
waveform.CPUsec = cputime - InitialCPUTime;
