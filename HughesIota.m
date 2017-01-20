function hughes_iota_deg = HughesIota(iota_deg, a, e, p)
%
% hughes_iota_deg = HughesIota(iota_deg, a, e, p)
%
% Returns Scott Hughes' definition of inclination (iota_hughes) iven the geometrictly defined 
% orbit parameters. Used by another function to determined iota given iota_hughes, 
% along with (a, e, p).
%
% The default accuracy is 1e-12;
%
% Steve Drasco
% 27 September 2007
%

% Do not allow inputs outside of 0 < iota < 180 
if(iota_deg > 180)
  iota_deg = 180;
end
if(iota_deg < 0)
  iota_deg = 0;
end

% compute E L Q
[E, L, Q] = ELzQ(a, e, p, iota_deg);
hughes_iota_deg = acos( L / sqrt(L^2 + Q) ) * 180 / pi;
