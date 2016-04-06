function sn = SNlisaCC(f,kovertyear)
% SNlisaCC      Returns array of noise power spectral density
%
% Synopsis:     sn = SNlisaCC(f)
%
% Input:        Array of frequencies, kappa/T (year^-1, default 1.5)
%               kappa/T = 0 yields no GWDB noise
%
% Output:       Array of real numbers
%
% Notes:        From Curt's notes

oneyear = 365.25 * 86400;

if nargin < 2
    kovert = 1.5 / oneyear;
else
    kovert = kovertyear / oneyear;
end
    
af = abs(f);

sinst = 0.75 * (9.18e-52 * af.^(-4) + 1.59e-41 + 9.18e-38 * af.^2);

sgal = 2.1e-45 * af.^(-7/3);

dndf = 2e-3 * af.^(-11/3);

sn = min(sinst .* exp(kovert * dndf), sinst + sgal);

% no DC power
sn(1) = sn(2);
