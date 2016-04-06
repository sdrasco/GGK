% fiducial choice
a = 0.5;
e0 = 0.8;
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
tspan = 35.15;
SmallSteps = 1258291;
coords = 'spherical';
BigSteps = 3;
tol = 1e-3;
D = 1.e17;
thetasb_deg = 45;
phisb_deg = 150;
theta_k_deg = 60;
phi_k_deg = 200;
order = 'quadrupole';
phase0 = [0 0 0 0];

% 
S = phase_InitializeWaveform(a, e0, p0, iota0_deg, phase0, M, mu, tspan, SmallSteps, BigSteps, tol, coords);
robs = 1;
thetaobs = 45;
phiobs = 0;
[hplus hcross]=ObserveWaveform(S, robs, thetaobs, phiobs, order);