function [hplus hcross]=hphx(hin, theta, phi, order)
%
% [hplus hcross]=hphx(hin, theta, phi, order)
%
% Computes h_plus and h_cross from cartesian components of 
% the radiation field h and the sky-osition of the observer.
%
% Input is:
%      hin = output structure from CartesianMetric.m
%    theta = polar angle of observers position (in degrees)
%      phi = azimuthal angle of observers position (in degrees)
%    order = string telling to use quardrupole approximation ('quadrupole')
%            or quadrupole plus octupole ('octupole')
%
% Output is:
%
%    hpluss = h_+ * r/mu
%    hcross = h_{\times} * r/mu
%
% See also CARTESIANMETRIC 
%
% Steve Drasco
%

% convert to radians
theta = theta * pi/180;
phi = phi * pi/180;

% make sure that order has been set
if ~strcmp(order,'quadrupole') && ~strcmp(order,'octupole')
	error('Order input must be set to either ''quadrupole'' or ''octupole''.');
end

if strcmp(order,'quadrupole')
    % use only the quadrupole formula
    h.xx = hin.quad.xx;
    h.yy = hin.quad.yy;
    h.zz = hin.quad.zz;
    h.xy = hin.quad.xy;
    h.xz = hin.quad.xz;
    h.yz = hin.quad.yz;

elseif strcmp(order,'octupole')
    % compute unit vector from source to observer;
    nx = 2*sin(theta) * cos(phi);
    ny = 2*sin(theta) * sin(phi);
    nz = 2*cos(theta);

    % project the octupole term, and sum with quadrupole term
    h.xx = hin.quad.xx + nx.*hin.oct.xxx + ny.*hin.oct.yxx + nz.*hin.oct.zxx;
    h.yy = hin.quad.yy + nx.*hin.oct.xyy + ny.*hin.oct.yyy + nz.*hin.oct.zyy;
    h.zz = hin.quad.zz + nx.*hin.oct.xzz + ny.*hin.oct.yzz + nz.*hin.oct.zzz;
    h.xy = hin.quad.xy + nx.*hin.oct.xxy + ny.*hin.oct.yxy + nz.*hin.oct.zxy;
    h.xz = hin.quad.xz + nx.*hin.oct.xxz + ny.*hin.oct.yxz + nz.*hin.oct.zxz;
    h.yz = hin.quad.yz + nx.*hin.oct.xyz + ny.*hin.oct.yyz + nz.*hin.oct.zyz;
end

% we'll use these a bunch of times, and they are a pain to retype each time
cts = cos(theta)^2;
sts = sin(theta)^2;
cps = cos(phi)^2;
sps = sin(phi)^2;

% theta-theta component
hTT = cts * (h.xx*cps + h.xy*sin(2*phi) + h.yy*sps ) ...
    + h.zz*sts - sin(2*theta)*(h.xz*cos(phi) + h.yz *sin(phi));

% theta-phi component
hTP = cos(theta)*(h.xy*cos(2*phi) + 0.5*h.yy*sin(2*phi) - 0.5*h.xx*sin(2*phi)) ...
    + sin(theta)*(h.xz*sin(phi) - h.yz*cos(phi));

% phi-phi component
hPP = h.xx*sps - h.xy*sin(2*phi) + h.yy*cps;

% plus and cross
hplus = 0.5*(hTT - hPP);
hcross = hTP;

% make them collumn vectors
hplus = hplus(:);
hcross = hcross(:);
