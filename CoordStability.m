

% load parameters
FishTest;

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

% record parameters
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
S.mu = mu ;
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
S.t0 = t0;

% take more BigSteps
S.BigSteps = S.BigSteps;
S.SmallSteps = SmallSteps;

% set Reomer delay in seconds if using LISA 
S.ReomerDelay_sec = 2000;
%S.ReomerDelay_sec = 0;

% get sky position of source as viewed from earth
[S.thetaE_deg S.phiE_deg] = EarthSkyPosition(S);

% set epsilon and delta
eps = 1e-4;
delta = 1e-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute derivatives at first point  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% store mu
output.a.mu = S.mu;

% build plus structure
Sp = S;
Sp.mu = S.mu*(1.0+eps);
Wp = InitializeWaveform(Sp);
      
% build minus structure
Sm = S;
Sm.mu = S.mu*(1.0-eps);
Wm = InitializeWaveform(Sm);

% just in case eps isn't well represented memory
TrueEps = 0.5*(Sp.mu - Sm.mu)/S.mu;

% get common time array for splines
tmin = max([min(Wp.x.t) min(Wm.x.t)]);
tmax = min([max(Wp.x.t) max(Wm.x.t)]);
t = linspace(tmin,tmax,length(Wm.x.t));

% compute derivatives
derivs1.dphidmu = (spline(Wp.x.t,Wp.x.phi,t) - spline(Wm.x.t,Wm.x.phi,t)) / (2 * TrueEps);
derivs1.ddphidmu = (spline(Wp.x.t,Wp.x.dphi,t) - spline(Wm.x.t,Wm.x.dphi,t)) / (2 * TrueEps);
derivs1.dddphidmu = (spline(Wp.x.t,Wp.x.ddphi,t) - spline(Wm.x.t,Wm.x.ddphi,t)) / (2 * TrueEps);
derivs1.drdmu = (spline(Wp.x.t,Wp.x.r,t) - spline(Wm.x.t,Wm.x.r,t)) / (2 * TrueEps);
derivs1.ddrdmu = (spline(Wp.x.t,Wp.x.dr,t) - spline(Wm.x.t,Wm.x.dr,t)) / (2 * TrueEps);
derivs1.dddrdmu = (spline(Wp.x.t,Wp.x.ddr,t) - spline(Wm.x.t,Wm.x.ddr,t)) / (2 * TrueEps);
derivs1.t = t;

% make sure nothing gets reused by mistake
clear Sp Sm Wp Wm t
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute derivatives at second point %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% step to second point
S.mu = S.mu*(1.0+delta);

% store mu
output.b.mu = S.mu;

% build plus structure
Sp = S;
Sp.mu = S.mu*(1.0+eps);
Wp = InitializeWaveform(Sp);
      
% build minus structure
Sm = S;
Sm.mu = S.mu*(1.0-eps);
Wm = InitializeWaveform(Sm);

% just in case eps isn't well represented memory
TrueEps = 0.5*(Sp.mu - Sm.mu)/S.mu;

% get common time array for splines
tmin = max([min(Wp.x.t) min(Wm.x.t)]);
tmax = min([max(Wp.x.t) max(Wm.x.t)]);
t = linspace(tmin,tmax,length(Wm.x.t));

% compute derivatives
derivs2.dphidmu = (spline(Wp.x.t,Wp.x.phi,t) - spline(Wm.x.t,Wm.x.phi,t)) / (2 * TrueEps);
derivs2.ddphidmu = (spline(Wp.x.t,Wp.x.dphi,t) - spline(Wm.x.t,Wm.x.dphi,t)) / (2 * TrueEps);
derivs2.dddphidmu = (spline(Wp.x.t,Wp.x.ddphi,t) - spline(Wm.x.t,Wm.x.ddphi,t)) / (2 * TrueEps);
derivs2.drdmu = (spline(Wp.x.t,Wp.x.r,t) - spline(Wm.x.t,Wm.x.r,t)) / (2 * TrueEps);
derivs2.ddrdmu = (spline(Wp.x.t,Wp.x.dr,t) - spline(Wm.x.t,Wm.x.dr,t)) / (2 * TrueEps);
derivs2.dddrdmu = (spline(Wp.x.t,Wp.x.ddr,t) - spline(Wm.x.t,Wm.x.ddr,t)) / (2 * TrueEps);
derivs2.t = t;

% make sure nothing gets reused by mistake
clear Sp Sm Wp Wm t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spline all derivatives onto common times %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time
tmin = max([min(derivs1.t) min(derivs2.t)]);
tmax = min([max(derivs1.t) max(derivs2.t)]);
output.t = linspace(tmin,tmax,length(derivs1.t));

% first point
output.a.dphidmu = spline(derivs1.t,derivs1.dphidmu,output.t);
output.a.ddphidmu = spline(derivs1.t,derivs1.ddphidmu,output.t);
output.a.dddphidmu = spline(derivs1.t,derivs1.dddphidmu,output.t);
output.a.drdmu = spline(derivs1.t,derivs1.drdmu,output.t);
output.a.ddrdmu = spline(derivs1.t,derivs1.ddrdmu,output.t);
output.a.dddrdmu = spline(derivs1.t,derivs1.dddrdmu,output.t);


% second point
output.b.dphidmu = spline(derivs2.t,derivs2.dphidmu,output.t);
output.b.ddphidmu = spline(derivs2.t,derivs2.ddphidmu,output.t);
output.b.dddphidmu = spline(derivs2.t,derivs2.dddphidmu,output.t);
output.b.drdmu = spline(derivs2.t,derivs2.drdmu,output.t);
output.b.ddrdmu = spline(derivs2.t,derivs2.ddrdmu,output.t);
output.b.dddrdmu = spline(derivs2.t,derivs2.dddrdmu,output.t);





