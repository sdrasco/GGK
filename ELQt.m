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
r = rp_ra(E,L,Q);
ra = r(1);
rp = r(2);
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