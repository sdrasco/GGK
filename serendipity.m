% script to try and reproduce the long ~180 cycle NR waveform described in Bela et al PRL (arxiv:1502.04953)

% initial parameters
a = 0; % dimensionless spin of big hole
e0 =7.5*34e-5; % long-run NR paper uses 34e-5.  too small for us.
p0 = 27.87966;  % gives M Omega_Phi_0. matching NR paper value of 6.7930e-3. 
iota0_deg = 15; % initial orbital inclination
theta0_deg = 90; % initial polar position
phi0_deg = 0; % initial azimuthal angle 
sign_rdot0 = 1; % initial sign of dr/dt
sign_Tdot0 = -1; % initial sign of dtheta/dt
m1 = 45.5/(1 + 1/7);  % Large mass in solar masses 
m2 = m1/7; % small mass in solar masses
M = m1 + m2; % total mass in solar masses (NR paper uses 45.5)
mu = m1*m2/M; % reduced mass
tspan = 30;  % this is in seconds.  NR paper shows ~ 23-24 seconds.
SmallSteps = 1e5; % number of time steps for waveform
coords = 'spherical'; % choice for how to treat Boyer-Lindquist coords
BigSteps = false;  % false to ignore, number if you want that many geodesics
tol = 1e-6; % accuracy used in KerrGeodesic (orbital parameters, spectra)
phase0 = [1 2 3 4];  % initial orbital phases, no reason for this choice

% main calculation here (computes each geodesic and stitches them together)
S = phase_InitializeWaveform(a, e0, p0, iota0_deg, phase0, M, mu, tspan, SmallSteps, BigSteps, tol, coords);

% clear everything that is now in S
clear a e0 p0 iota0_deg r0 theta0_deg phi0_deg sign_rdot0 sign_Tdot0 ...
    M mu tspan SmallSteps coords BigSteps tol order phase0

% observe the waveform
h.t = S.x.t;
h.robs = 1;
h.thetaobs = 45; % theta location of observer on hole's sky
h.phiobs = 0; % phi location of observer on hole's sky
h.order = 'quadrupole';
[h.plus h.cross]=ObserveWaveform(S, h.robs, h.thetaobs, h.phiobs, h.order);
