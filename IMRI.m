% common parameters
M = 1e6;
mu = 1e3;
tspan = 3.1556926e7;
SmallSteps = 2^20;
BigSteps = 2^5;
tol = 1e-4;
coords = 'spherical';
D = 1.e17;
thetasb_deg = 45;
phisb_deg = 150;
theta_k_deg = 60;
phi_k_deg = 200;
order = 'octupole';
t0=0;
phi0_deg = 29;
sign_rdot0 = -1;
sign_Tdot0 = -1;

% spin 0.1 runs
% a = 0.1;
% e0 = 4.9903e-02;
% p0 =2.4317e+01;
% iota0_deg = 4.4905e+01;
% r0 = 0.99*p0/(1-e0);
% theta0_deg = 1.01*(90 - iota0_deg);
% S = FisherMatrix(a, e0, p0, iota0_deg, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D);
% save 'S_0.1_0.01.mat' S
% clear S a e0 p0 iota0_deg r0 theta0_deg;
% %
% a = 0.1;
% e0 = 5.6232e-01;
% p0 =2.3433e+01;
% iota0_deg = 4.4895e+01;
% r0 = 0.99*p0/(1-e0);
% theta0_deg = 1.01*(90 - iota0_deg);
% S = FisherMatrix(a, e0, p0, iota0_deg, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D);
% save 'S_0.1_0.1.mat' S
% clear S a e0 p0 iota0_deg r0 theta0_deg;
%
%%% this block gives a problem
% a = 0.1;
% e0 = 9.1377e-01;
% p0 = 1.8887e+01;
% iota0_deg = 4.4880e+01;
% r0 = 0.99*p0/(1-e0);
% theta0_deg = 1.01*(90 - iota0_deg);
% S = FisherMatrix(a, e0, p0, iota0_deg, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D);
% save 'S_0.1_0.25.mat' S
% clear S a e0 p0 iota0_deg r0 theta0_deg;


% spin 0.5 runs
% a = 0.5;
% e0 = 8.5689e-02;   
% p0 = 2.3967e+01;
% iota0_deg = 4.4375e+01;
% r0 = 0.99*p0/(1-e0);
% theta0_deg = 1.01*(90 - iota0_deg);
% S = FisherMatrix(a, e0, p0, iota0_deg, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D);
% save 'S_0.5_0.01.mat' S
% clear S a e0 p0 iota0_deg r0 theta0_deg;
% %
% a = 0.5;
% e0 = 7.6804e-01;
% p0 = 2.1488e+01;
% iota0_deg = 4.4257e+01;
% r0 = 0.99*p0/(1-e0);
% theta0_deg = 1.01*(90 - iota0_deg);
% S = FisherMatrix(a, e0, p0, iota0_deg, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D);
% save 'S_0.5_0.1.mat' S
% clear S a e0 p0 iota0_deg r0 theta0_deg;
%

%%% this block gives a problem
% a = 0.5;
% e0 = 9.9796e-01;
% p0 =1.0663e+01;   
% iota0_deg = 4.4307e+01;
% r0 = 0.99*p0/(1-e0);
% theta0_deg = 1.01*(90 - iota0_deg);
% S = FisherMatrix(a, e0, p0, iota0_deg, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D);
% save 'S_0.5_0.25.mat' S
% clear S a e0 p0 iota0_deg r0 theta0_deg;

% spin 0.8 runs
% a = 0.8;
% e0 = 1.7300e-01;
% p0 = 2.3712e+01;
% iota0_deg = 4.3806e+01;
% r0 = 0.99*p0/(1-e0);
% theta0_deg = 1.01*(90 - iota0_deg);
% S = FisherMatrix(a, e0, p0, iota0_deg, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D);
% save 'S_0.8_0.01.mat' S
% clear S a e0 p0 iota0_deg r0 theta0_deg;
%

a = 0.8;
e0 = 9.6874e-01;
p0 =1.5545e+01;
iota0_deg = 4.3213e+01;
r0 = 0.99*p0/(1-e0);
theta0_deg = 1.01*(90 - iota0_deg);
S = FisherMatrix(a, e0, p0, iota0_deg, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D);
save 'S_0.8_0.1.mat' S
clear S a e0 p0 iota0_deg r0 theta0_deg;
%
a = 0.8;
e0 = 9.9867e-01;
p0 = 9.7699e+00;
iota0_deg = 4.3137e+01;
r0 = 0.99*p0/(1-e0);
theta0_deg = 1.01*(90 - iota0_deg);
S = FisherMatrix(a, e0, p0, iota0_deg, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D);
save 'S_0.8_0.25.mat' S
clear S a e0 p0 iota0_deg r0 theta0_deg;


