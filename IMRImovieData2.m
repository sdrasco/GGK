% fiducial choice
a = 0.9;
e0 = 0.5;
p0 =6;
iota0_deg = 50;
r0 = p0 / ( 1 - e0);
theta0_deg = 90;
phi0_deg = 0;
sign_rdot0 = 1;
sign_Tdot0 = -1;
M = 3000;
mu = 1.4;
t0=0;
%tspan = 17.15;
tspan = 643.0725;
SmallSteps = 1258291;
coords = 'spherical';
BigSteps = 64;
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