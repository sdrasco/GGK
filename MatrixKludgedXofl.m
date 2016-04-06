function x=MatrixKludgedXofl(inspiral,lambda, lambda_x0)
%
% x=MatrixKludgedXofl(inspiral,lambda)
%
% An attempt at a quicker version of KludgedXofl.  This version splits the 
% calculation into chunks in an effort to minimize memory consumption, and uses
% matrices for some of the for loops.  Tests seem to show that, despite my efforts,
% this code is actually slower than a more simple minded code (KLUDGEDXOFL).
%
% Returns Boyer-Lindquist (BL) coordinates and pseudo-cartesian (PC) coordinates
%
%         x.t = BL t/M
%         x.r = BL r/M
%     x.theta = BL theta
%       x.phi = BL phi
%    x.CPUsec = number of seconds spent building this structure
%
% as a function of lambda, given the following inputs:
%
%  inspiral: output structure from KLUDGEDGEODESIC
%    lambda: array of lambda values at which you want positions
%
% Output also has first, second, and third derivatives of the world 
% line.  Note that t is differentiated with respect to lambda, and 
% while spatial coordinates are differentiated with respect to t.
%
%         x.dt = dt/dlambda
%        x.ddt = d/dlambda(dt/dlambda)
%       x.dddt = d/dlambda(d/dlambda(dt/dlambda))
%         x.dr = dr/dt
%        x.ddr = d/dt(dr/dt)
%       x.dddr = d/dt(d/dt(dr/dt))
%     x.dtheta = dtheta/dt
%    x.ddtheta = d/dt(dtheta/dt)
%   x.dddtheta = d/dt(d/dt(dtheta/dt))
%       x.dphi = dphi/dt
%      x.ddphi = d/dt(dphi/dt)
%     x.dddphi = d/dt(d/dt(dphi/dt))
%
% See also KLUDGEDXOFL KERRGEODESIC
%
% Steve Drasco
%

% start clock to keep track of computational cost
InitialCPUTime = cputime;

if nargin < 3
  % assume a fiducial initial phase convention
  %lambda_x0 = [0 0 inspiral.geodesic{1}.LT/2 0];
  lambda_x0 = [0 0 0 0];
end

% make sure that lambda is a column vector
x.lambda = lambda(:);

% initialize some memory to speed things
% since Matlab complains that loops take longer if 
% variables are growing in size as the loop proceeds
x.t = zeros(size(lambda));
x.dt = zeros(size(lambda));
x.ddt = zeros(size(lambda));
x.dddt = zeros(size(lambda));
x.r = zeros(size(lambda));
dr = zeros(size(lambda));
ddr = zeros(size(lambda));
dddr = zeros(size(lambda));
x.theta = zeros(size(lambda));
dtheta = zeros(size(lambda));
ddtheta = zeros(size(lambda));
dddtheta = zeros(size(lambda));
x.phi = zeros(size(lambda));
dphi = zeros(size(lambda));
ddphi = zeros(size(lambda));
dddphi = zeros(size(lambda));

% Memory costs are high if lambda gets too long.
% Split calculation into smaller blocks if necessary.  
% Speed tests on JPL desktop (2.5 GHz dual G5) home laptop (1.5 GHz G4) 
% did best when MaxLength was about 300-400 (home-JPL).  
MaxLength = 400;
lblocks = ceil(length(lambda) / MaxLength);
for lindex = 1:lblocks
  if lindex < lblocks
    lrange = ((lindex-1)*MaxLength + 1):lindex*MaxLength;
  else
    lrange = ((lindex-1)*MaxLength + 1):length(x.lambda);
  end
  lambda = x.lambda(lrange);

  % evolve positional elements by requiring dw/dlambda = Upsilon
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
  clear Ek;

  % build sine and cosine matrices for n-expansions
  En = []; %En = zeros(length(lambda),nmax);
  for n=1:nmax;
    En = [En Er.^n]; % En(:,n) = Er.^n;
  end
  cn = real(En);
  sn = -imag(En);
  clear En;

  % evaluate splines to build matrix of coordinate expansion coefficients
  tn=zeros(nmax,length(lambda));
  dtn=zeros(nmax,length(lambda));
  ddtn=zeros(nmax,length(lambda));
  dddtn=zeros(nmax,length(lambda));
  rn=zeros(nmax,length(lambda));
  drn=zeros(nmax,length(lambda));
  ddrn=zeros(nmax,length(lambda));
  dddrn=zeros(nmax,length(lambda));
  phin=zeros(nmax,length(lambda));
  dphin=zeros(nmax,length(lambda));
  ddphin=zeros(nmax,length(lambda));
  dddphin=zeros(nmax,length(lambda));
  for n=1:nmax
    tn(n,:) = ppval(inspiral.tn_PP{n},lambda);
    dtn(n,:) = n*tn(n,:);
    ddtn(n,:) = n*dtn(n,:);
    dddtn(n,:) = n*ddtn(n,:);
    rn(n,:) = ppval(inspiral.rn_PP{n+1},lambda); % inspiral.rn starts from n=0
    drn(n,:) = n*rn(n,:);
    ddrn(n,:) = n*drn(n,:);
    dddrn(n,:) = n*ddrn(n,:);
    phin(n,:) = ppval(inspiral.phin_PP{n},lambda);
    dphin(n,:) = n*phin(n,:);
    ddphin(n,:) = n*dphin(n,:);
    dddphin(n,:) = n*ddphin(n,:);
  end
  tk=zeros(kmax,length(lambda));
  dtk=zeros(kmax,length(lambda));
  ddtk=zeros(kmax,length(lambda));
  dddtk=zeros(kmax,length(lambda));  
  phik=zeros(kmax,length(lambda));
  dphik=zeros(kmax,length(lambda));
  ddphik=zeros(kmax,length(lambda));
  dddphik=zeros(kmax,length(lambda));
  Tk=zeros(kmax,length(lambda));
  dTk=zeros(kmax,length(lambda));
  ddTk=zeros(kmax,length(lambda));
  dddTk=zeros(kmax,length(lambda));
  for k=1:kmax
    tk(k,:) = ppval(inspiral.tk_PP{k},lambda);
    dtk(k,:) = k*tk(k,:);
    ddtk(k,:) = k*dtk(k,:);
    dddtk(k,:) = k*ddtk(k,:);
    phik(k,:) = ppval(inspiral.phik_PP{k},lambda);
    dphik(k,:) = k*phik(k,:);
    ddphik(k,:) = k*dphik(k,:);
    dddphik(k,:) = k*ddphik(k,:);
    Tk(k,:) = ppval(inspiral.thetak_PP{k+1},lambda); % inspiral.Tk starts from n=0
    dTk(k,:) = k*Tk(k,:);
    ddTk(k,:) = k*dTk(k,:);
    dddTk(k,:) = k*ddTk(k,:);
  end

  % evaluate splines for r_0 and theta_0
  r0=ppval(inspiral.rn_PP{1},lambda); 
  T0=ppval(inspiral.thetak_PP{1},lambda);

  % compute BL coordinates and first three derivatives
  x.t(lrange) = wt + diag(sn*tn) + diag(sk*tk);
  x.dt(lrange) = ppval(inspiral.Gamma_PP,lambda) ...
      + ppval(inspiral.UpsilonR_PP,lambda).*diag(cn*dtn) ...
      + ppval(inspiral.UpsilonTheta_PP,lambda).*diag(ck*dtk); 
  x.ddt(lrange) = -(ppval(inspiral.UpsilonR_PP,lambda).^2).*diag(sn*ddtn) ...
      - (ppval(inspiral.UpsilonTheta_PP,lambda).^2).*diag(sk*ddtk);
  x.dddt(lrange) = -(ppval(inspiral.UpsilonR_PP,lambda).^3).*diag(cn*dddtn) ...
      - (ppval(inspiral.UpsilonTheta_PP,lambda).^3).*diag(ck*dddtk);
  x.r(lrange) = r0 + 2*diag(cn*rn); 
  dr(lrange) = -2*ppval(inspiral.UpsilonR_PP,lambda).*diag(sn*drn); 
  ddr(lrange) = -2*(ppval(inspiral.UpsilonR_PP,lambda).^2).*diag(cn*ddrn); 
  dddr(lrange) = 2*(ppval(inspiral.UpsilonR_PP,lambda).^3).*diag(sn*dddrn); 
  x.theta(lrange) = T0 + 2*diag(ck*Tk);
  dtheta(lrange) = -2*ppval(inspiral.UpsilonTheta_PP,lambda).*diag(sk*dTk); 
  ddtheta(lrange) = -2*(ppval(inspiral.UpsilonTheta_PP,lambda).^2).*diag(ck*ddTk); 
  dddtheta(lrange) = 2*(ppval(inspiral.UpsilonTheta_PP,lambda).^3).*diag(sk*dddTk); 
  x.phi(lrange) = wp + diag(sn*phin) + diag(sk*phik);
  dphi(lrange) = ppval(inspiral.UpsilonPhi_PP,lambda) ...
      + ppval(inspiral.UpsilonR_PP,lambda).*diag(cn*dphin) ...
      + ppval(inspiral.UpsilonTheta_PP,lambda).*diag(ck*dphik); 
  ddphi(lrange) = -(ppval(inspiral.UpsilonR_PP,lambda).^2).*diag(sn*ddphin) ...
      - (ppval(inspiral.UpsilonTheta_PP,lambda).^2).*diag(sk*ddphik);
  dddphi(lrange) = -(ppval(inspiral.UpsilonR_PP,lambda).^3).*diag(cn*dddphin) ...
      - (ppval(inspiral.UpsilonTheta_PP,lambda).^3).*diag(ck*dddphik);
   
  % probably not necessary, but clear some things by hand in case matlab
  % does a bad job with memory management.
  clear ck sk cn sn tn dtn ddtn dddtn rn drn ddrn dddrn phin dphin ...
      ddphin dddphin phik dphik ddphik dddphik Tk dTk ddTk dddTk;

end

% The derivatives computed above are lambda-derivatives.  For the spatial
% coodinates, convert them to t-derivatives now. We will need first five 
% inverse powers of dt repeatedly. 
dtm1 = x.dt.^(-1);
dtm2 = x.dt.^(-2);
dtm3 = x.dt.^(-3);
dtm4 = x.dt.^(-4);
dtm5 = x.dt.^(-5);

% lambda-derivatives to t-derivatives, freeing space as we go
x.dr = dr .* dtm1;
x.ddr = (ddr .* dtm2) - (dr .* x.ddt .* dtm3);
x.dddr = (dddr .* dtm3) - (3 * ddr .* x.ddt .* dtm4) ...
    - (dr .* x.dddt .* dtm4) + (3 * dr .* x.ddt.^2 .* dtm5);
clear dr ddr dddr;
x.dtheta = dtheta .* dtm1;
x.ddtheta = (ddtheta .* dtm2) - (dtheta .* x.ddt .* dtm3);
x.dddtheta = (dddtheta .* dtm3) - (3 * ddtheta .* x.ddt .* dtm4) ...
    - (dtheta .* x.dddt .* dtm4) + (3 * dtheta .* x.ddt.^2 .* dtm5);
clear dtheta ddtheta dddtheta;
x.dphi = dphi .* dtm1;
x.ddphi = (ddphi .* dtm2) - (dphi .* x.ddt .* dtm3);
x.dddphi = (dddphi .* dtm3) - (3 * ddphi .* x.ddt .* dtm4) ...
    - (dphi .* x.dddt .* dtm4) + (3 * dphi .* x.ddt.^2 .* dtm5);
clear dphi ddphi dddphi;

% Clear things just in case Matlab is lazy about it, since they're big
clear dtm1 dtm2 dtm3 dtm4

% log the computational cost of this job
x.CPUsec = cputime - InitialCPUTime;
