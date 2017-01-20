function [thetaE_deg phiE_deg] = EarthSkyPosition(S)
%
% [thetaE_deg phiE_deg] = EarthSkyPosition(S)
%
% compute thetaE and phiE: direction to Earth in frame of BH, where BH 
% spin is around its z-axis.  phiE is azimuthal location of Earth wrt 
% MBH spin.  We have freedom to set this equal to zero, which we do for 
% simplicity.  Then phi0 is phi coord of particle at t0, WITH RESPECT TO 
% DIRECTION TO EARTH -- i.e., the diff in these 2 angles.
%

% check and read inputs
theta_k = S.theta_k_deg * pi/180;
phi_k = S.phi_k_deg * pi/180;
thetasb = S.thetasb_deg * pi/180;
phisb = S.phisb_deg * pi/180;

% main calculation
kdotn = cos(theta_k)*cos(thetasb) + sin(theta_k)*sin(thetasb)*cos(phisb-phi_k);
thetaE = acos(-kdotn);
phiE = 0;
thetaE_deg = thetaE*180.0/pi;
phiE_deg = phiE*180.0/pi;