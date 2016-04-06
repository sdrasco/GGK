function [hI hII t] = LISAResponse(W, S, hplus, hcross, RoemerDelay_sec)
%
% Adapted from Curt's LISAh2
%
% Steve Drasco (May, 2009)


% note: these first two blocks are not just to check the structure, but to 
% remind me what parameters this function is actually using.

% check input structure W 
if ~isfield(W,'SecPerM') 
    error('LISAResponse: ''SecPerM'' field missing from input structure W');
end
if ~isfield(W,'x') 
    error('LISAResponse: ''x'' field missing from input structure W');
end
if ~isfield(W.x,'t') 
    error('LISAResponse: ''t'' field missing from input structure W.x');
end
if ~isfield(W,'t0') 
    error('LISAResponse: ''t0'' field missing from input structure W');
end



% check input structure S 
if ~isfield(S,'tspan') 
    error('LISAResponse: ''tspan'' field missing from input structure S');
end
if ~isfield(S,'SmallSteps') 
    error('LISAResponse: ''SmallSteps'' field missing from input structure S');
end
if ~isfield(S,'thetasb_deg') 
    error('LISAResponse: ''thetasb_deg'' field missing from input structure S');
end
if ~isfield(S,'phisb_deg') 
    error('LISAResponse: ''phisb_deg'' field missing from input structure S');
end
if ~isfield(S,'theta_k_deg') 
    error('LISAResponse: ''theta_k_deg'' field missing from input structure S');
end
if ~isfield(S,'phi_k_deg') 
    error('LISAResponse: ''phi_k_deg'' field missing from input structure S');
end
if ~isfield(S,'ReomerDelay_sec')
    error('LISAResponse: ''ReomerDelay_sec'' field is missing from the input structure S');
end

% get long t-vector in seconds
t = W.SecPerM * W.x.t;

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

thetasb = S.thetasb_deg * pi/180;
phisb = S.phisb_deg * pi/180;
theta_k = S.theta_k_deg * pi/180;
phi_k = S.phi_k_deg * pi/180;
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
hplus_shift = spline(t,hplus,t_shift);
hcross_shift = spline(t,hcross,t_shift);

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

% Resample data with even spacing.  tt covers the range of t (with 
% RoemerDelay clipped off both ends) but is evenly sampled
Wt0sec = W.t0 * W.SecPerM;
%tI = min([Wt0sec Wt0sec+S.tspan]);
%tF = max([Wt0sec Wt0sec+S.tspan]);
tI = min([Wt0sec Wt0sec+S.tspan]) + S.ReomerDelay_sec;
tF = max([Wt0sec Wt0sec+S.tspan]) - S.ReomerDelay_sec;
tt = linspace(tI,tF,S.SmallSteps);
%tt = linspace(min(t_shift),max(t_shift),S.SmallSteps);
hI = spline(t_shift, hIu, tt);
hII = spline(t_shift, hIIu, tt);

% output is called t, not tt
t = tt;

end