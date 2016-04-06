function S = STest(eps,tol,fname)
% Tests the code in mfiles--code generating numerical kludge waveforms,
% their inner prods, Fisher matrices, statistical and model error, etc.

% start tracking CPU cost of construction
InitialCPUtime = cputime;

%here are some fiducial param values:
a = 0.6;
e0 = 0.3;
%p0 =10;
p0=4.894670474209899e+00;
iota_deg0 = 30;
%r0 = p0*(1 +e0/3.0); 
r0 = p0 / (1+e0);
theta0_deg = 70;
phi0_deg = 112;
sign_rdot0 = 1;
sign_Tdot0 = -1;
M = 1.e6;
mu = 10.e0;


%short waveforms for quick testing
%tI = -1.e5; t0 = 0.e0; tF = 1.e5; 
%SmallSteps = 1024*16;
%BigSteps = 8;
%tol = 0.1;

%long waveforms:
%tI = -1.5e7; t0 = 0.e0; tF = 1.e7;
tI = -(365.25*2*3600); t0 = 0; tF = 0;
SmallSteps = 1024*1024;
BigSteps = 64;
%tol = 1e-4;

delt = (tF -tI)/(SmallSteps -1);
coords = 'spherical';

%quasi-random distance and angles, for testing
%D = 1.e17;
%thetasb = 45*pi/180.0;
%phisb = 150*pi/180.0;
%theta_k = 60*pi/180.0;
%phi_k = 200*pi/180.0;

% mimic source from Barack and Cutler's table III
D = (1.e17) * (54/30) * (34.1/30);
thetasb = pi/4;
phisb = 0;
theta_k = pi/8;
phi_k = 0;



%Uncomment this block to calculate the snr-squared
%[hI hII] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
%       sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
%hIf = cfft(hI,tI,delt,SmallSteps);
%hIIf = cfft(hII,tI,delt,SmallSteps);
%snr2= innerprod(hIf,hIf,hIIf,hIIf,delt,SmallSteps);

nparam = 14;
[DhIf DhIIf] = SDerivs(eps,nparam, a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
        sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
%Now calculating Fisher matrix
%Coords are lnM, ln mu, a(dimensionless spin), p0 (semilatus rectum in units of M), e0, iota_deg0, 3 phases, ln r, theta, phi (11 in all)

fishmat = zeros(nparam,nparam);
%Now making Fisher matrix
for i = 1:nparam
    for j = 1:i
        display(['Computing F_{' num2str(i) ',' num2str(j) '} ...']);
        fishmat(i,j) = innerprod(DhIf(i,:),DhIf(j,:),DhIIf(i,:),DhIIf(j,:),delt,SmallSteps);
        fishmat(j,i) = fishmat(i,j);
    end
end
%L = chol(fishmat, 'lower');  %This is correct version of call for Matlab
%version on Curt's desktop
L = chol(fishmat);  % Curt's laptop's version of matlab chol function takes a single function
invL = inv(L);
%invfishC = transpose(inv(L))*inv(L);    % <-------------------------------
invfishC = invL*transpose(invL);
% parameter derivatives are wrt (1-14):
%1: ln D
%2: ln M
%3: ln mu
%4: a \equiv S/M^2
%5: p0 \equiv p/M at time t0
%6  r0 \equiv r/M at time t0
%7: e0
%8: iota0 in radians
%9: theta_0 in radians
%10: phi0 (radians)
%11: thetasb in radians
%12:phisb (radians)
%13: theta_k in radians
%14:phi_k (radians)
S.dlogD = sqrt(invfishC(1,1));
S.dlogM = sqrt(invfishC(2,2));
S.dlogmu = sqrt(invfishC(3,3));
S.da = sqrt(invfishC(4,4));
S.dp0 = sqrt(invfishC(5,5));
S.dr0 = sqrt(invfishC(6,6));
S.de0 = sqrt(invfishC(7,7));
S.diota0 = sqrt(invfishC(8,8));
S.dtheta0 = sqrt(invfishC(9,9));
S.dphi0 = sqrt(invfishC(10,10));
S.dthetasb = sqrt(invfishC(11,11));
S.dphisb = sqrt(invfishC(12,12));
S.dtheta_k = sqrt(invfishC(13,13));
S.dphi_k = sqrt(invfishC(14,14));
S.delomega = 2.e0*pi*sqrt(invfishC(11,11)*invfishC(12,12) - invfishC(11,12)*invfishC(11,12));
S.a = a;
S.e0 = e0;
S.p0 = p0;
S.iota_deg0 = iota_deg0;
S.r0 = r0;
S.theta0_deg = theta0_deg;
S.phi0_deg = phi0_deg;
S.sign_rdot0 = sign_rdot0;
S.sign_Tdot0 = sign_Tdot0;
S.M = M;
S.mu = mu;
S.tol = tol;
S.eps = eps;
disp(['dlogD=', num2str(S.dlogD)]);
disp(['dlogM=', num2str(S.dlogM)]);
disp(['dlogmu=', num2str(S.dlogmu)]);
disp(['da=', num2str(S.da)]);
disp(['dp0=', num2str(S.dp0)]);
disp(['dphi0=', num2str(S.dphi0)]);
disp(['dthetasb=', num2str(S.dthetasb)]);
disp(['dphisb=', num2str(S.dphisb)]);
disp(['dtheta_k=', num2str(S.dtheta_k)]);
disp(['dphi_k=', num2str(S.dphi_k)]);
disp(['delomega(sky error box in steradians)=', num2str(S.delomega)]);


% send results to a file
S.F = fishmat;
S.iF = invfishC;
S.CPUsec = cputime - InitialCPUtime;
eval(['save ''' fname ''' S']);


disp('stopping here for debugging purposes');
disp('stopped');
