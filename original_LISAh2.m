function [hI hII] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
        sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k)

% implement Cutler (low-f) response for LISA channel I and II
% based on "LISApatterns.m" from cleanlisabias directory

%like LISAh.m. except in this version:
%t0 is some fiducial time when the waveform sweeps past some frequency at the SS Barycenter
%tI and tII are time-endpts at the LISA detector

% - theta, phi: sky position of the source (ecliptic), in radians

% - phib (LISA's angular position around Sun) = phib0 at t=0 (not
% necessarily at t = t0).
%   Below we set phib0 to be 0, for simplicity. 

%we make waveforms going foward and backward from t0.
%AU = 499.00478370;
%solar_mass = (AU^3.0)*(0.01720209895/86400)^2.0;
%using Drasco's less accurate version for now, for consistency
SecPerMsun = 4.9255e-6;

BigStepsI = floor(BigSteps*(t0-tI)/(tF-tI));
BigStepsF = BigSteps - BigStepsI;
SmallStepsI = floor(SmallSteps*(t0-tI)/(tF-tI)); 
SmallStepsF = SmallSteps - SmallStepsI;


%Now we make grid from which we'll later interpolate a bit denser
SmallStepsI = floor(SmallStepsI*1.2);
SmallStepsF = floor(SmallStepsF*1.2);

%Rem we later have to account for fact that both I and F
%solutions give h(t0), so we have a "duplicate" ppoint

%We extend a little
%beyond tI and tF in each direction so we can account for th Reomer delay,
%which is a maximum of 500 sec. (We extend by 600 just to be safe.) 
%extra = 2000;
extra = 0; % TEST

% Steve's InitialzeWaveform says it wants t0 (5th argument) in units of M
t0_M = t0/(M*SecPerMsun);

waveformI = OriginalInitializeWaveform(a, e0, p0, iota_deg0, t0_M, r0, theta0_deg, phi0_deg, ...
              sign_rdot0, sign_Tdot0, M, mu, (tI-extra -t0), SmallStepsI, BigStepsI, tol, coords);
%waveformI = OriginalInitializeWaveform(a, e0, p0, 180.0-iota_deg0, t0_M, r0, 180.0-theta0_deg, phi0_deg + 180.0, ...
%              -1.0*sign_rdot0, -1.0*sign_Tdot0, M, mu, abs(tI-extra -t0), SmallStepsI, BigStepsI, tol, coords);
waveformF = OriginalInitializeWaveform(a, e0, p0, iota_deg0, t0_M, r0, theta0_deg, phi0_deg, ... 
              sign_rdot0, sign_Tdot0, M, mu, tF+extra -t0, SmallStepsF, BigStepsF, tol, coords);

%Next construct hplus and hcross in the barycenter

%thetaE and phiE are direction to Earth in frame of BH, where BH spin is
%around its z-axis:

kdotn = cos(theta_k)*cos(thetasb) + sin(theta_k)*sin(thetasb)*cos(phisb-phi_k);
thetaE = acos(-kdotn);
%phiE is azimuthal location of Earth wrt MBH spin.  We have freedom to set
%this equal to zero, which we do for simplicity.   Then phio is phi coord
%of particle at t0, WITH RESPECT TO DIRECTION TO EARTH -- i.e., the diff in
%these 2 angles.
phiE = 0;
%rem Steve routine wants theta and phi
% in degrees.  We also convert phiE, in case we ever want to change
% the convention that its zero
thetaE_deg = thetaE*180.0/pi;
phiE_deg = phiE*180.0/pi;

[hplusI hcrossI]=ObserveWaveform(waveformI, D, thetaE_deg, phiE_deg, 'quadrupole');
[hplusF hcrossF]=ObserveWaveform(waveformF, D, thetaE_deg, phiE_deg, 'quadrupole');
%[hplusI hcrossI]=ObserveWaveform(waveformI, D, thetaE_deg, phiE_deg, 'octupole');
%[hplusF hcrossF]=ObserveWaveform(waveformF, D, thetaE_deg, phiE_deg, 'octupole');


hplusI = hplusI';
hcrossI = hcrossI';
hplusF = hplusF';
hcrossF = hcrossF';

%putting earlies times at the far left:
hplusI = fliplr(hplusI);
hplusI = hplusI(1:SmallStepsI-1);
%cutting off the midpoint of whole integration, since it is in both the
%I(initial) and F(final) parts
hcrossI = fliplr(hcrossI);
hcrossI = hcrossI(1:SmallStepsI-1);
tIvec = waveformI.x.t;
tIvec = fliplr(tIvec);
tIvec = tIvec(1:SmallStepsI-1);
tFvec = waveformF.x.t;
t = waveformI.SecPerM*[tIvec,tFvec];
hplus = [hplusI,hplusF];
hcross = [hcrossI,hcrossF];

%Steve's ObserveWaveform gives (r/mu)*h; below corrects for that
%hplus = (mu*SecPerMsun/D)*hplus;
%hcross =(mu*SecPerMsun/D)*hcross;

%saving some memory:
clear waveformI; clear waveformF;
clear hplusI; clear hplusF; clear hcrossI; clear hcrossF;
clear tIvec; clear tFvec;

% LISA c.o.m. trajectory
% initial azimuth of LISA
alpha0 = 0;
phib0 = 0;

years = 365.25 * 86400; %is this approximation consistent with other formulae in FI, FII?
Omega = 2 * pi / years;
%as test, setting Omega =0 and phib0 =1.0
%Omega = 0.e0;
%phib0 = 1.0;


%phib = phib0 + Omega * (t-t0);
phib = phib0 + Omega * t;  %we set phib (location of LISA) to be phib0 at t=0 (NOT necessarily at t = t0)

% LISA angles

%                             this 1 is i
%alpha1 = Omega * (t -t0) - pi/12 - (1 - 1)*(pi/3) + alpha0;
alpha1 = Omega *t  - pi/12 - (1 - 1)*(pi/3) + alpha0;  %(specifying LISA orientation at t =0, 
%which is not necess at t = t0).

%thetas and phis are source position in LISA frame, while thetasb and phisb
%are source position in ecliptic frame (with "b" at end standing for
%"barycenter").

costhetas = 0.5 * cos(thetasb) - (sqrt(3)/2) * sin(thetasb) * cos(phib - phisb);
phis = alpha1 + pi/12 + atan2( sqrt(3)*cos(thetasb) + sin(thetasb)*cos(phib - phisb), ...
                               2*sin(thetasb)*sin(phib - phisb) );                           
        
%Now calculating polarization angle psi--see Cutler notes.
%Calculation based on fact that Gair convention for hplus and hcross has x
%along theta direction and y along phi, in spherical coords centered on
%MBH, with MBH spin along noth pole

%First calculating p and q directions (see notes):
%Again, k is direction of MBH spin in ecliptic coords, and n is direction
%from Earth to source

nvec = [sin(thetasb)*cos(phisb), sin(thetasb)*sin(phisb), cos(thetasb)];
kvec = [sin(theta_k)*cos(phi_k), sin(theta_k)*sin(phi_k), cos(theta_k)];
kdotn = dot(kvec,nvec);

CON = (1 - kdotn^2.0)^0.5;
qvec = cross(nvec,kvec)/CON;
pvec = (-kvec + (kdotn)*nvec)/CON;

%zL(t) is normal to LISA plane
zLz = 0.5*ones(1,length(t));
zLx = -0.5*sqrt(3)*cos(phib);
zLy = -0.5*sqrt(3)*sin(phib);

zLdotq = qvec(1)*zLx + qvec(2)*zLy + qvec(3)*zLz;
zLdotp = pvec(1)*zLx + pvec(2)*zLy + pvec(3)*zLz;
%zLdotp 

psis = atan2(zLdotq, zLdotp);

FpI  = 0.5 * (1 + costhetas.^2) .* cos(2*phis) .* cos(2*psis) ...
                   - 1.0*costhetas .* sin(2*phis) .* sin(2*psis);

%for checking: 
%FpIch  = 0.5 * (1 + costhetas.^2) .* cos(2*phis) .* cos(2*psis); 
%FpIch = FpIch - costhetas .* sin(2*phis) .* sin(2*psis);
               
FxI  = 0.5 * (1 + costhetas.^2) .* cos(2*phis) .* sin(2*psis) ...
                   + costhetas .* sin(2*phis) .* cos(2*psis);


% For II, using the phis -> phis - pi/4 transformation [ cos(2x) ->  sin(2x),
%                                                sin(2x) -> -cos(2x) ]
% certainly this stands to be optimized (or is MATLAB already doing it?)...

FpII = 0.5 * (1 + costhetas.^2) .* sin(2*phis) .* cos(2*psis) ...
                   + costhetas .* cos(2*phis) .* sin(2*psis);

               
FxII = 0.5 * (1 + costhetas.^2) .* sin(2*phis) .* sin(2*psis) ...
                   - costhetas .* cos(2*phis) .* cos(2*psis);

% Now we shift hplus and hcross by Roemer delay:
%Use this:
% Doppler phase (2 * pi * RoverC = 3135.35 s)

%Now we multiply hplus and hcross by LISA response functions, making use of fact that 
%Drasco is using Gair&Glampedakis convention for hplus and hcross, which has opposite handedness
%of Cutler convention; in practice Cutler has to multiply Drasco hcross by -1
%to get back to Cutler convention
AU = 149597870691/299792458;

%checking setting AU =0;
%AU = 0;  -doesn't affect disagreement with fig in Cutler 98.

t_shift = t - AU*sin(thetasb) .* cos(phib - phisb);

%now cutting 600 sec off left and right ends of t_shift, so can interpolate onto it
%from t-grid
delt1 = t(2) - t(1);
dN1 = floor(600/delt1);
sizet = length(t);
delt2 = t(sizet) - t(sizet-1);
dN2 = floor(600/delt2);
%now cutting off a little from left and right ends  
t_shift = t_shift(dN1+1:sizet-dN2);

%now interpolating to find values of h_plus and h_cross at times t_shift
hplus_shift = interp1(t,hplus,t_shift,'spline');
%hcross_shift = interp1(t,hcross,t_shift,'spline','extrap')
hcross_shift = interp1(t,hcross,t_shift,'spline');

%Now see what happens if we just switch sign of hcross
%hcross_shift = -hcross_shift;
%hcross_shift = 0;

%tplot = t(dN1+1:SmallSteps-dN2);
%figure;
%plot(tplot,hplus,t,hplus_shift);
%hold on

%Fp and Fx are unshifted, since they're naturally written in terms of
%Barycenter time
FpI = FpI(dN1+1:sizet-dN2);
FxI = FxI(dN1+1:sizet-dN2);
FpII = FpII(dN1+1:sizet-dN2);
FxII = FxII(dN1+1:sizet-dN2);
hIu  = (sqrt(3)/2)*(FpI.*hplus_shift + FxI.*hcross_shift);
hIIu  = (sqrt(3)/2)*(FpII.*hplus_shift + FxII.*hcross_shift);

%saving some memory:
clear FpI; clear FxI; clear FpII; clear FxII;

%now sampliong data with even spacing
%tt covers same range as t, but is evenly sampled

%tt = 0:SmallSteps-1;
%tt = tI + tt*(tF-tI)/(SmallSteps-1);
tt = linspace(min(t_shift),max(t_shift),SmallSteps);

hI = spline(t_shift, hIu, tt);
%figure;
%plot(tt,hI);
hII = spline(t_shift, hIIu, tt);
%figure;
%plot(tt,hII);
%saving some memory:
clear hIu; clear hIIu; 

%disp('stopping here for debugging purposes');
%disp('stopped');
