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
%above from my unpublished "timing.pdf" 
t0 = t0/solar_mass;
tf = tf/solar_mass;

%options = odeset('RelTol', 1e-4); % value from curt's older uses
options = odeset('RelTol', 1e-6);
tspan = linspace(t0,tf,Ntimes);

% solve once with equally spaced time values
[trange, solution] = ode45(@ELQdots,tspan,[E0 L0 Q0], options);


% determine a better spacing in time (to get evenly spaced energy values)
%
% Curt had this in here. Not sure needed -- doubling the work?
% We could spline to get a better spacing, but if we had a very course
% spacing to begin with, that could introduce significant error.
%
Et = solution(:,1);
Erange = linspace(Et(1), Et(end), Ntimes);
tspan = spline(Et, trange, Erange);

% solve again with the improved spacing in time
[trange, solution] = ode45(@ELQdots,tspan,[E0 L0 Q0], options);
Et = solution(:,1);
Lt = solution(:,2);
Qt = solution(:,3);

%converting trange from solar mass to sec
trange = trange*solar_mass;

%figure
%plot(trange,Et/m)
%figure
%plot(trange,Lt/(m*M))
%figure
%plot(trange,Qt/(m*m*M*M))
%disp('stopping here');
%disp('stopped');


