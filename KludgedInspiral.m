function inspiral = KludgedInspiral(a, M, mu, t, e, p, iota_deg, t0, tol)
%
% inspiral = KludgedInspiral(a, M, mu, t, e, p, iota_deg, t0, tol)
%
% Given the dimensionless parameters
%
%         t = dimensionless coordinate time (t/M)
%         a = dimensionless black hole spin parameter ( 0 <= a < 1)
%         e = eccentricity 
%         p = dimensionless semilatus-rectum
%  iota_deg = inclination, in degrees [defined as: pi/2 - sgn(L) theta_min]
%        t0 = t(lambda=0) (M)
%       tol = requested fractional accuracy
%
% this program produces a structure "inspiral" that contains 
% a collection of data that can describes other details of 
% the orbit.
%
% Steve Drasco
% 30 May 2008
%

% start tracking CPU cost of construction
InitialCPUTime = cputime;

% initialize empty matrices and vectors that will be indexed later
tn=[];
tk=[];
rn=[];
phik=[];
phin=[];
Tk=[];
Gamma=zeros(size(t));
Ur=zeros(size(t));
UT=zeros(size(t));
Up=zeros(size(t));

% copy main parameters into inspiral structure
inspiral.a = a;
inspiral.M = M;
inspiral.mu = mu;
SecPerMsun = 4.9255e-6;
inspiral.SecPerM = SecPerMsun * inspiral.M;
inspiral.SecPermu = SecPerMsun * inspiral.mu;
inspiral.lambda = zeros(size(t));
inspiral.t = t;
inspiral.tol = tol;

% build an orbit structure for each set of orbital constants
for i=1:length(t)

    % compute geodesic structure
    inspiral.geodesic{i} = KerrGeodesic(a, e(i), p(i), iota_deg(i), tol);

    % status report on time.
    if i == 1
        SecSpent = [];
    end
    SecSpent = [SecSpent; inspiral.geodesic{i}.CPUsec];
    SecLeft =  round((length(t) - i) * mean(SecSpent));
    display(['Computed Geodesic ' num2str(i) '/' num2str(length(t)) ...
        '.  About ' num2str(SecLeft) ' sec remaining.']);
    
    % build matrices of coefficients to be used for spline fits
    tn = [tn inspiral.geodesic{i}.tn];
    tk = [tk inspiral.geodesic{i}.tk];
    rn = [rn inspiral.geodesic{i}.rn];
    Tk = [Tk inspiral.geodesic{i}.thetak];
    phin = [phin inspiral.geodesic{i}.phin];
    phik = [phik inspiral.geodesic{i}.phik];
    Gamma(i)=inspiral.geodesic{i}.Gamma;
    Ur(i)=inspiral.geodesic{i}.UpsilonR;
    UT(i)=inspiral.geodesic{i}.UpsilonTheta;
    Up(i)=inspiral.geodesic{i}.UpsilonPhi;    
end

% compute lambda at each value of t
% (sampled on the long time scale)
ooG = 1./Gamma;
ooG_pp = spline(inspiral.t,ooG);
IooG_pp = fnint(ooG_pp);
inspiral.lambda = ppval(IooG_pp,inspiral.t) - ppval(IooG_pp,t0);

% compute spline polynomials for the frequencies
inspiral.Gamma_PP = spline(inspiral.lambda,Gamma);
inspiral.UpsilonR_PP = spline(inspiral.lambda,Ur);
inspiral.UpsilonTheta_PP = spline(inspiral.lambda,UT);
inspiral.UpsilonPhi_PP = spline(inspiral.lambda,Up);

% compute spine polynomials for integrals of frequencies
inspiral.IYt_pp = fnint(inspiral.Gamma_PP);
inspiral.IYr_pp = fnint(inspiral.UpsilonR_PP);
inspiral.IYT_pp = fnint(inspiral.UpsilonTheta_PP);
inspiral.IYp_pp = fnint(inspiral.UpsilonPhi_PP);

% get index maxima
nmax = inspiral.geodesic{1}.nmax;
kmax = inspiral.geodesic{1}.kmax;

% compute spline polynomials for the series expansion coefficients
for n=1:nmax
    inspiral.tn_PP{n} = spline(inspiral.lambda,tn(n,:));
    inspiral.rn_PP{n} = spline(inspiral.lambda,rn(n,:));
    inspiral.phin_PP{n} = spline(inspiral.lambda,phin(n,:));
end
for k=1:kmax
    inspiral.tk_PP{k} = spline(inspiral.lambda,tk(k,:));
    inspiral.thetak_PP{k} = spline(inspiral.lambda,Tk(k,:));
    inspiral.phik_PP{k} = spline(inspiral.lambda,phik(k,:));
end

% theta and r expansions started with index values of zero.  Do spline fit
% for their last elements "by hand"
inspiral.rn_PP{nmax+1} = spline(inspiral.lambda,rn(nmax+1,:));
inspiral.thetak_PP{kmax+1} = spline(inspiral.lambda,Tk(kmax+1,:));

% record computational cost
inspiral.CPUsec = cputime - InitialCPUTime;
