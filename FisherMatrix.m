function S = FisherMatrix(a, e0, p0, iota0_deg, t0, ...
    r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, ...
    M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, ...
    thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D, ParameterKey)
%
% S = FisherMatrix(a, e0, p0, iota0_deg, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D)
%
% Or if you want to use a limited set of parameters, use 
%
% S = FisherMatrix(a, e0, p0, iota0_deg, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0, M, mu, tspan, SmallSteps, BigSteps, tol, coords, order, thetasb_deg, phisb_deg, theta_k_deg, phi_k_deg, D, ParameterKey)
%
% Computes the Fisher Matrix for an EMRI. 
%
% Takes the following inputs:
%
%    [to be filled in later]
%
% The ParameterKey input is optional.  If unspecified, the default is used:
%
%    ParameterKey={
%        'ln D', ...
%        'ln M', ...
%        'ln mu', ...
%        'a', ...
%        'p0', ...
%        'r0', ...
%        'e0', ...
%        'iota0', ...
%        'theta0', ...
%        'phi0', ...
%        'thetasb', ...
%        'phisb', ...
%        'theta_k', ...
%        'phi_k', ...
%    };
%
% If ParameterKey is specified, the ordering of its elements are not 
% important, but each of its elments must be one of those in the default 
% set.
% 
% Output structure has the following fields:
%
%    [to be filled in later]
% 
% Adapted from Curt's test routine, test.m.
%
% Steve Drasco (May 2009)
%

% start tracking CPU cost of construction
InitialCPUtime = cputime;

% sanity checks for inputs
if log2(SmallSteps) ~= round(log2(SmallSteps))
    original_SmallSteps = SmallSteps;
    SmallSteps = 2^round(log2(SmallSteps));
    disp(['FisherMatrix: We require SmallSteps to be a power of ' ...
        '2. Rounding SmallSteps from ' num2str(original_SmallSteps) ...
        ' to ' num2str(SmallSteps) '.']);
end
if log2(BigSteps) ~= round(log2(BigSteps))
    original_BigSteps = BigSteps;
    BigSteps = 2^round(log2(BigSteps));
    disp(['FisherMatrix: We require BigSteps to be a power of ' ...
        '2. Rounding BigSteps from ' num2str(original_BigSteps) ...
        ' to ' num2str(BigSteps) '.']);
end

% Define parameter key if it has not been entered as an input
if nargin == 23
    ParameterKey={
        'ln D', ...
        'ln M', ...
        'ln mu', ...
        'a', ...
        'p0', ...
        'r0', ...
        'e0', ...
        'iota0', ...
        'theta0', ...
        'phi0', ...
        'thetasb', ...
        'phisb', ...
        'theta_k', ...
        'phi_k', ...
    };
end

% verify that the parameter key obeys the naming rules
nparam = length(ParameterKey);
for i=1:nparam
    if strcmp(ParameterKey{i},'ln D')
    elseif strcmp(ParameterKey{i},'ln M')
    elseif strcmp(ParameterKey{i},'ln mu')
    elseif strcmp(ParameterKey{i},'a')
    elseif strcmp(ParameterKey{i},'p0')
    elseif strcmp(ParameterKey{i},'r0')
    elseif strcmp(ParameterKey{i},'e0')
    elseif strcmp(ParameterKey{i},'iota0')
    elseif strcmp(ParameterKey{i},'theta0')
    elseif strcmp(ParameterKey{i},'phi0')
    elseif strcmp(ParameterKey{i},'thetasb')
    elseif strcmp(ParameterKey{i},'phisb')
    elseif strcmp(ParameterKey{i},'theta_k')
    elseif strcmp(ParameterKey{i},'phi_k')
    else
        error(['FisherMatrix: unrecognoied ParameterKey ' ...
            ParameterKey{i}]);
    end
end


% record input arguments
S.a = a;
S.e0 = e0;
S.p0 = p0;
S.iota0_deg = iota0_deg;
S.r0 = r0;
S.theta0_deg = theta0_deg;
S.phi0_deg = phi0_deg;
S.sign_rdot0 = sign_rdot0;
S.sign_Tdot0 = sign_Tdot0;
S.M = M;
S.mu = mu;
S.tspan = tspan;
S.SmallSteps = SmallSteps;
S.BigSteps = BigSteps;
S.tol = tol;
S.coords = coords;
S.order = order;
S.thetasb_deg = thetasb_deg;
S.phisb_deg = phisb_deg;
S.theta_k_deg = theta_k_deg;
S.phi_k_deg = phi_k_deg;
S.D = D;
S.ParameterKey = ParameterKey;
S.nparam = length(ParameterKey);
S.delt = abs(tspan) / (S.SmallSteps-1);
S.t0 = t0;

% set Reomer delay in seconds if using LISA 
S.ReomerDelay_sec = 2000;
%S.ReomerDelay_sec = 0;

% get sky position of source as viewed from earth
[S.thetaE_deg S.phiE_deg] = EarthSkyPosition(S);

% compute derivatives
[DhIf DhIIf] = GenericDerivs(S);

% calculate Fisher matrix
S.Gamma = zeros(nparam,nparam);
for i = 1:nparam
    for j = 1:i
        S.Gamma(i,j) = innerprod(DhIf(i,:),DhIf(j,:),...
                                 DhIIf(i,:),DhIIf(j,:),S.delt,S.SmallSteps);
        S.Gamma(j,i) = S.Gamma(i,j);
    end
end

% invert Fisher matrix
L = chol(S.Gamma);  
invL = inv(L);
S.invGamma = invL*transpose(invL);

% Record error estimates as indicidual fields. This block is written in
% such a way allow for incomplete or jumbled parameter keys.
for i=1:nparam
    if strcmp(ParameterKey{i},'ln D')
        S.dlogD = sqrt(S.invGamma(i,i));
        S.SNR = sqrt(S.Gamma(i,i));
    elseif strcmp(ParameterKey{i},'ln M')
        S.dlnM = sqrt(S.invGamma(i,i));
    elseif strcmp(ParameterKey{i},'ln mu')
        S.dlnmu = sqrt(S.invGamma(i,i));
    elseif strcmp(ParameterKey{i},'a')
        S.da = sqrt(S.invGamma(i,i));
    elseif strcmp(ParameterKey{i},'p0')
        S.dp0 = sqrt(S.invGamma(i,i));
    elseif strcmp(ParameterKey{i},'r0')
        S.dr0 = sqrt(S.invGamma(i,i));
    elseif strcmp(ParameterKey{i},'e0')
        S.de0 = sqrt(S.invGamma(i,i));
    elseif strcmp(ParameterKey{i},'iota0')
        S.diota0 = sqrt(S.invGamma(i,i));
    elseif strcmp(ParameterKey{i},'theta0')
        S.dtheta0 = sqrt(S.invGamma(i,i));
    elseif strcmp(ParameterKey{i},'phi0')
        S.dphi0 = sqrt(S.invGamma(i,i));
    elseif strcmp(ParameterKey{i},'thetasb')
        S.dthetasb = sqrt(S.invGamma(i,i));
        tsbIndex = i;
    elseif strcmp(ParameterKey{i},'phisb')
        S.dphisb = sqrt(S.invGamma(i,i));
        psbIndex = i;
    elseif strcmp(ParameterKey{i},'theta_k')
        S.dtheta_k = sqrt(S.invGamma(i,i));
    elseif strcmp(ParameterKey{i},'phi_k')
        S.dphi_k = sqrt(S.invGamma(i,i));
    else
        error(['FisherMatrix: unrecognoied ParameterKey:' ...
            ParameterKey{i}]);
    end
end
if isfield(S,'dthetasb') && isfield(S,'dphisb')
    S.delomega = 2*pi ...
     *sqrt(S.invGamma(tsbIndex,tsbIndex)*S.invGamma(psbIndex,psbIndex) ...
     - S.invGamma(tsbIndex,psbIndex)*S.invGamma(tsbIndex,psbIndex));
end

% record total CPU time
S.CPUsec = cputime - InitialCPUtime;

end