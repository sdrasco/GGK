function p = innermost(a, e, iota_deg)
%
% p = innermost(a, e, iota_deg)
%
% Returns the smallest value of p such that an orbit characterized by (a,
% e, iota_deg, p) will be stable [positive value of radial derivative of 
% radial potential].
%
% Steve Drasco (May, 2009)
%

% set accuracy request
options = optimset('TolX',1e-6);

pGuess = 3;
p = fzero(@(p) dVrdr(p, a, e, iota_deg), pGuess, options);

end

function output = dVrdr(p, a, e, iota_deg)
% output = dVrdr(p, a, e, iota_deg)
%
%   Returns derivative of radial potential with respect to r.  Semilatus
%   rectum (p) can be a vector input.
%
% Steve Drasco (May, 2009)

% simple parameters
rmin = p ./ (1+e);
varpi2 = rmin.^2 + a.^2;
Delta = varpi2 - 2*rmin;

% harder parameters
E = zeros(size(p));
L = zeros(size(p));
Q = zeros(size(p));
for i=1:length(p);
    [E, L, Q] = ELzQ(a, e, p(i), iota_deg);
end

% compute d[V_r(rmin)]/dr
output = 4 * E .* rmin .* (E .* varpi2 - a .* L);
output = output - (2*rmin - 2) .* ( rmin.^2 + (L-a.*E).^2 + Q );
output = output - 2*rmin.*Delta;

end


