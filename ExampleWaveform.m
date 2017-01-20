% fiducial choice
a = 0.5;
e0 = 0.1;
p0 =10;
iota0_deg = 33;
r0 = p0 / ( 1 - e0);
theta0_deg = 90;
phi0_deg = 0;
sign_rdot0 = 1;
sign_Tdot0 = -1;
M = 250;
mu = 1.4;
t0=0;
tspan = 59.71;
SmallSteps = 1e5;
coords = 'spherical';
BigSteps = 20;
tol = 1e-9;
D = 1.e17;
thetasb_deg = 45;
phisb_deg = 150;
theta_k_deg = 60;
phi_k_deg = 200;
phase0 = [0 0 0 0];

% main calculation here (computes each geodesic and stitches them together)
S = phase_InitializeWaveform(a, e0, p0, iota0_deg, phase0, M, mu, tspan, SmallSteps, BigSteps, tol, coords);

% clear everything that is now in S
clear a e0 p0 iota0_deg r0 theta0_deg phi0_deg sign_rdot0 sign_Tdot0 M mu t0 tspan SmallSteps coords BigSteps tol D thetasb_deg phisb_deg theta_k_deg phi_k_deg order phase0

% observe the waveform
h.t = S.x.t;
h.robs = 1;
h.thetaobs = 45;
h.phiobs = 0;
h.order = 'quadrupole';
[h.plus h.cross]=ObserveWaveform(S, h.robs, h.thetaobs, h.phiobs, h.order);