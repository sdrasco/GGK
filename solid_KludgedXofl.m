function x=solid_KludgedXofl(inspiral,lambda)
%
% x=solid_KludgedXofl(inspiral,lambda)
%
% Returns Boyer-Lindquist (BL) coordinates and pseudo-cartesian (PC) coordinates
%
%         x.t = BL t/M
%         x.r = BL r/M
%     x.theta = BL theta
%       x.phi = BL phi
%         x.x = PC x/M = x.r sin(x.theta) cos(x.phi)
%         x.y = PC y/M = x.r sin(x.theta) sin(x.phi)
%         x.z = PC z/M = x.r cos(x.theta)
%    x.lambda = Mino time * M = lambda
%    x.CPUsec = number of seconds spent building this structure
%
% as a function of lambda, given the following inputs:
%
%  inspiral: output structure from KLUDGEDGEODESIC
%    lambda: array of lambda values at which you want positions
%
%
% See also KLUDGEDGEODESIC KERRGEODESIC
%
% Steve Drasco
% 31 May 2008
%

% start clock to keep track of computational cost
InitialCPUTime = cputime;

% assume a fiducial initial phase convention
lambda_x0 = [0 0 inspiral.geodesic{1}.LT/2 0];

% make sure that lambda is a column vector
lambda = lambda(:);

% compute w's in one of two ways
%
% way #1: require that lambda_x is constant
%x.Gamma = ppval(inspiral.Gamma_PP,lambda);
%x.UpsilonR = ppval(inspiral.UpsilonR_PP,lambda);
%x.UpsilonTheta = ppval(inspiral.UpsilonTheta_PP,lambda);
%x.UpsilonPhi = ppval(inspiral.UpsilonPhi_PP,lambda);
%wt =        x.Gamma.*(lambda - lambda_x0(1));  
%wr =     x.UpsilonR.*(lambda - lambda_x0(2));  
%wT = x.UpsilonTheta.*(lambda - lambda_x0(3));  
%wp =   x.UpsilonPhi.*(lambda - lambda_x0(4));  
%
%
% way #2: require dw/dlambda = Upsilon
wt0 = -inspiral.geodesic{1}.Gamma * lambda_x0(1);
wr0 = -inspiral.geodesic{1}.UpsilonR * lambda_x0(2);
wT0 = -inspiral.geodesic{1}.UpsilonTheta * lambda_x0(3);
wp0 = -inspiral.geodesic{1}.UpsilonPhi * lambda_x0(4);
wt = wt0 + ppval(inspiral.IYt_pp, lambda);
wr = wr0 + ppval(inspiral.IYr_pp, lambda);
wT = wT0 + ppval(inspiral.IYT_pp, lambda);
wp = wp0 + ppval(inspiral.IYp_pp, lambda);

% compute exponentials
Er = exp(-i*wr);
ET = exp(-i*wT);

% get index maxima for series
nmax = inspiral.geodesic{1}.nmax;
kmax = inspiral.geodesic{1}.kmax;

% build sine and cosine matrices for k-expansions
Ek = []; % Ek = zeros(length(lambda),kmax);
for k=1:kmax;
  Ek = [Ek ET.^k]; % Ek(:,k) = ET.^k;
end
ck = real(Ek);
sk = -imag(Ek);

% build sine and cosine matrices for n-expansions
En = []; %En = zeros(length(lambda),nmax);
for n=1:nmax;
  En = [En Er.^n]; % En(:,n) = Er.^n;
end
cn = real(En);
sn = -imag(En);


% evaluate splines to build matrix of coordinate expansion coefficients
tn=zeros(nmax,length(lambda));
rn=zeros(nmax,length(lambda));
phin=zeros(nmax,length(lambda));
for n=1:nmax
    tn(n,:) = ppval(inspiral.tn_PP{n},lambda);
    rn(n,:) = ppval(inspiral.rn_PP{n+1},lambda); % inspira.rn starts from n=0
    phin(n,:) = ppval(inspiral.phin_PP{n},lambda);
end
tk=zeros(kmax,length(lambda));
phik=zeros(kmax,length(lambda));
Tk=zeros(kmax,length(lambda));
for k=1:kmax
    tk(k,:) = ppval(inspiral.tk_PP{k},lambda);
    phik(k,:) = ppval(inspiral.phik_PP{k},lambda);
    Tk(k,:) = ppval(inspiral.thetak_PP{k+1},lambda); % inspira.Tk starts from n=0
end

% evaluate splines for r_0 and theta_0
r0=ppval(inspiral.rn_PP{1},lambda); 
T0=ppval(inspiral.thetak_PP{1},lambda);

% compute BL coordinates
x.t = wt + diag(sn*tn) + diag(sk*tk);
x.r = r0 + 2*diag(cn*rn); 
x.theta = T0 + 2*diag(ck*Tk); 
x.phi = wp + diag(sn*phin) + diag(sk*phik);

% compute SC coordinates
x.x = x.r .* sin(x.theta) .* cos(x.phi);  % note: thought BL coords generalize to spheroidal
x.y = x.r .* sin(x.theta) .* sin(x.phi);  % rather than spherical coordinates.  need to 
x.z = x.r .* cos(x.theta);                % double check this.
x.lambda = lambda;
x.CPUsec = cputime - InitialCPUTime;

