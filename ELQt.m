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

% set accuracy tollerance (older code used 1e-4 without disaster)
options = odeset('RelTol', 1e-6,'Events',@StopCondition);

% solve the ODEs
sol = ode45(@ELQdots,[t0 tf],[E0 L0 Q0], options);
trange = sol.x;
Et = sol.y(1,:);
Lt = sol.y(2,:);
Qt = sol.y(3,:);

% If wanted, restrict to specified number of geodesics (otherwise Ntimes
% aka BigSteps will be ignored)
% if length(sol.x) ~= Ntimes
%     tspan = linspace(sol.x(1),sol.x(end),Ntimes);
%     trange = tspan;
%     solution = deval(sol,trange); % uses ode45's native interpolation
%     Et = solution(1,:);
%     Lt = solution(2,:);
%     Qt = solution(3,:);
% else
%     trange = sol.x;
%     Et = sol.y(1,:);
%     Lt = sol.y(2,:);
%     Qt = sol.y(3,:);
% end

%converting trange from solar mass to sec
trange = trange*solar_mass;

end

function [value,isterminal,direction] = StopCondition(t,y)

global M;

E = y(1);
L = y(2); 
Q = y(3);
r = rp_ra(E,L,Q);
ra = r(1);
rp = r(2);
p = 2*ra*rp/(ra+rp);
e = (ra-rp)/(ra+rp);

% Would be great to have a good way to stopping the integration.  For now
% just set to some minimum p and or e.  If a better stopping
% condition isnt' determined, turn p_stop and e_stop into arguments and
% contol their values elsewhere.
p_stop = 6.6;
e_stop = 1.375e-2;
if p/M < p_stop || e < e_stop
   value = 0; % we met the condition, will stop integration now
else
   value = 1; % any number other than 0 to say we haven't met condition
end
isterminal = 1; % stop the integration (0 keeps going, records meetin condition)
direction = 0; % direction of approach to stopping condition (we don't use)
end