function [hplus hcross]=ObserveWaveform(waveform, r, theta, phi, order)
%
% [hplus hcross]=ObserveWaveform(waveform, r, theta, phi, order)
%
% Given location of observer in spherical coordinates
%
%      r: radial distance in units of M
%  theta: polar angle in degrees
%    phi: azimuthal angle in degrees
%  order: string telling to use quardrupole approximation ('quadrupole')
%         or quadrupole plus octupole ('octupole')
%
% as well as a waveform structure produced by InitializeWaveform, 
% this routine computes h_+ and h_x evaluated at the coordinate time or 
% lambda values given in
%
%    waveform.x.t
%    waveform.x.lambda
%
% respectively.
%
% NOTE: currently this routine returns the normalized waveform
%
%     (r/mu) [h_+ h_x]
%
% If you want the actual waveform, scale by mu/r.
%
% See also INITIALIZEWAVEFORM
%
% By: Steve Drasco
%

% Feed the waveform's Cartesian metric components into the hphx routine
[hplus hcross]=hphx(waveform.h, theta, phi, order);

% normalize according to distance and mass
hplus = hplus * (waveform.mu * waveform.SecPerMsun / r);
hcross = hcross *  (waveform.mu * waveform.SecPerMsun / r);

