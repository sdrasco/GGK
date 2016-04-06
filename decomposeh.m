function hmkn = decomposeh(a, e, p, iota_deg, Nperoid, coords, order, tol)
%
% just a test routine for now
%

% check inputs for errors
if ~strcmp(coords,'spherical') && ~strcmp(coords,'spheroidal')
	error('coords must be either ''spherical'' or ''spheroidal''');
end
if ~strcmp(order,'quadrupole') && ~strcmp(order,'octupole')
    error('order must be either ''quadrupole'' or ''octupole''');
end

% start tracking CPU cost of construction
InitialCPUtime = cputime;

% compute the orbit
orbit = KerrGeodesic(a, e, p, iota_deg, tol);
LongestPeriod = max([orbit.LP orbit.LT orbit.LR]);
lambda = linspace(1,Nperoid*LongestPeriod,1e4); 
r0 = orbit.rn(1) + 2*sum(orbit.rn(2:end)); 
theta0_deg= (orbit.thetak(1) + 2*sum(orbit.thetak(2:end))) * 180/pi; 
phi0_deg=0; 
t0=0; 
sign_rdot0=1; 
sign_Tdot0=1;
orbit.x = xofl(orbit, lambda, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0);

% get t-range for integrals
tmin = orbit.x.t(100);
tmax = orbit.x.t(end-100);
T = tmax - tmin;

% calculate cartesian components of the worldline and derivatives
if strcmp(coords,'spherical')
    % approximating Boyer-Lindquist as sphericlal flat coordinates.
    Z = CartesianWorldline_Spherical(orbit.x);
else
    % approximating Boyer-Lindquist as spheroidal flat coordinates.
    Z = CartesianWorldline_Spheroidal(orbit.x, orbit.a);
end
orbit.x.x = Z.x;
orbit.x.y = Z.y;
orbit.x.z = Z.z;
orbit.x.dx = Z.dx;
orbit.x.dy = Z.dy;
orbit.x.dz = Z.dz;
orbit.x.ddx = Z.ddx;
orbit.x.ddy = Z.ddy;
orbit.x.ddz = Z.ddz;
orbit.x.dddx = Z.dddx;
orbit.x.dddy = Z.dddy;
orbit.x.dddz = Z.dddz;

% compute cartesian components of the metric
orbit.h = CartesianMetric(orbit.x);

% compute waveform
theta = 37;
phi = 89;
[hplus hcross]=hphx(orbit.h, theta, phi, order);

% make a spline-version of h over integration range
hplus_PP = spline(orbit.x.t,hplus);
hcross_PP = spline(orbit.x.t,hcross);
splineh = @(t) ppval(hplus_PP,t) - (1i) * ppval(hcross_PP,t);

% loop to collect frequencies
orbit.mmax = max([orbit.kmax orbit.nmax]);
Nfreq = 5^3; %orbit.kmax * orbit.nmax * orbit.mmax;
FreqMatrix = zeros(Nfreq, 4);
ind = 0;
for m=1:5 %orbit.mmax
    for k=1:5 %orbit.kmax
        for n=1:5 % orbit.nmax
            ind = ind + 1;
            FreqMatrix(ind,:)=[m, k, n, ...
                m*orbit.OmegaPhi + k*orbit.OmegaTheta + n*orbit.OmegaR];
        end
    end
end

% options for integration
options = odeset('RelTol',1e-2);

% main loop to decompose h into hmkn
Y = eye(Nfreq);
I = zeros(Nfreq,1);
omega = FreqMatrix(:,4);
for i=1:Nfreq
    integrand = @(t,y) splineh(t) .* exp((1i)*omega(i)*t); % second argument needed for ode45
    [tt,II] = ode45(integrand,[tmin tmax], 0, options); I(i) = II(end)/T;
    %I(i) = quadl(integrand, tmin, tmax, 1e-2) / T;
    display([num2str(FreqMatrix(i,1:3)) ' : integral number ' num2str(i) '/' num2str(Nfreq) ' : ' num2str(I(i))]);
    for j=1:Nfreq
        if i~=j
            Y(i,j) = ( (1i)/(T*(omega(i) - omega(j))) ) ...
                * ( 1 - exp((1i)*(omega(i)-omega(j))*T) );
        end
    end
end

% invert mode-mixing-matrix to get the modes
modes = Y\I;

% store all sorts of things in the output structure
hmkn.FreqMatrix = FreqMatrix;
hmkn.modes = modes;
hmkn.T = T;
hmkn.Y = Y;
hmkn.I = I;
hmkn.orbit = orbit;

% record CPU time
hmkn.CPUsec = cputime - InitialCPUtime;