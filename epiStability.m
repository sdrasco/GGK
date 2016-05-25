% load parameters
FishTest;

% force number of steps to be powers of 2 (not needed in this test code, 
% but keeping it for consistancy)
if log2(SmallSteps) ~= round(log2(SmallSteps))
    SmallSteps = 2^round(log2(SmallSteps));
end
if log2(BigSteps) ~= round(log2(BigSteps))
    BigSteps = 2^round(log2(BigSteps));
end

% convert to hughes' iota, and store unit conversion parameters
iota_hughes_deg0 = HughesIota(iota0_deg, a, e0, p0);
SecPerMsun = 4.9255e-6;
SecPerM = SecPerMsun * M;

% set epsilon, used in derivative calculations
eps = 1e-6;

% set delta, used to tell where we evaluate derivatives
delta = 1e-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute derivatives at first point  (called A) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% store mu
output.A.mu = mu;

% plus parameters
mup = mu*(1.0+eps);
[tp, pp, ep, ip, Ep, Lp, Qp] = EvolveOrbit(p0,e0,iota_hughes_deg0*pi/180,...
   t0*SecPerM,t0*SecPerM + tspan,BigSteps,M,a*M,mup);
      
% build minus parameters
mum = mu*(1.0-eps);
[tm, pm, em, im, Em, Lm, Qm] = EvolveOrbit(p0,e0,iota_hughes_deg0*pi/180,...
    t0*SecPerM,t0*SecPerM + tspan,BigSteps,M,a*M,mum);

% just in case eps isn't well represented in memory
TrueEps = 0.5*(mup - mum)/mu;

% get common time array for splines
tmin = max([min(tp) min(tm)]);
tmax = min([max(tp) max(tm)]);
t = linspace(tmin,tmax,length(tp));

% compute derivatives
derivsA.dedmu = (spline(tp,ep,t) - spline(tm,em,t)) / (2 * TrueEps);
derivsA.dpdmu = (spline(tp,pp,t) - spline(tm,pm,t)) / (2 * TrueEps);
derivsA.didmu = (spline(tp,ip,t) - spline(tm,im,t)) / (2 * TrueEps);
derivsA.dEdmu = (spline(tp,Ep,t) - spline(tm,Em,t)) / (2 * TrueEps);
derivsA.dLdmu = (spline(tp,Lp,t) - spline(tm,Lm,t)) / (2 * TrueEps);
derivsA.dQdmu = (spline(tp,Qp,t) - spline(tm,Qm,t)) / (2 * TrueEps);

% store time in derivs structure too
derivsA.t = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute derivatives at second point  (called B) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% just be extra sure that nothing is left from last block of code
clear tp pp ep ip Ep Lp Qp tm pm em im Em Lm Qm 
clear mup mum TrueEps t tmin tmax

% step to second point
mu = mu*(1.0+delta);

% store mu
output.B.mu = mu;

% plus parameters
mup = mu*(1.0+eps);
[tp, pp, ep, ip, Ep, Lp, Qp] = peitELQ(p0,e0,iota_hughes_deg0*pi/180,...
    t0*SecPerM,t0*SecPerM + tspan,BigSteps,M,a*M,mup);
      
% build minus parameters
mum = mu*(1.0-eps);
[tm, pm, em, im, Em, Lm, Qm] = peitELQ(p0,e0,iota_hughes_deg0*pi/180,...
    t0*SecPerM,t0*SecPerM + tspan,BigSteps,M,a*M,mum);

% just in case eps isn't well represented in memory
TrueEps = 0.5*(mup - mum)/mu;

% get common time array for splines
tmin = max([min(tp) min(tm)]);
tmax = min([max(tp) max(tm)]);
t = linspace(tmin,tmax,length(tp));

% compute derivatives
derivsB.dedmu = (spline(tp,ep,t) - spline(tm,em,t)) / (2 * TrueEps);
derivsB.dpdmu = (spline(tp,pp,t) - spline(tm,pm,t)) / (2 * TrueEps);
derivsB.didmu = (spline(tp,ip,t) - spline(tm,im,t)) / (2 * TrueEps);
derivsB.dEdmu = (spline(tp,Ep,t) - spline(tm,Em,t)) / (2 * TrueEps);
derivsB.dLdmu = (spline(tp,Lp,t) - spline(tm,Lm,t)) / (2 * TrueEps);
derivsB.dQdmu = (spline(tp,Qp,t) - spline(tm,Qm,t)) / (2 * TrueEps);

% store time in derivs structure too
derivsB.t = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spline all derivatives onto common times %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time
tmin = max([min(derivsA.t) min(derivsB.t)]);
tmax = min([max(derivsA.t) max(derivsB.t)]);
output.t = linspace(tmin,tmax,length(derivsA.t));

% first point
output.A.dedmu = spline(derivsA.t,derivsA.dedmu,output.t);
output.A.dpdmu = spline(derivsA.t,derivsA.dpdmu,output.t);
output.A.didmu = spline(derivsA.t,derivsA.didmu,output.t);
output.A.dEdmu = spline(derivsA.t,derivsA.dEdmu,output.t);
output.A.dLdmu = spline(derivsA.t,derivsA.dLdmu,output.t);
output.A.dQdmu = spline(derivsA.t,derivsA.dQdmu,output.t);

% second point
output.B.dedmu = spline(derivsB.t,derivsB.dedmu,output.t);
output.B.dpdmu = spline(derivsB.t,derivsB.dpdmu,output.t);
output.B.didmu = spline(derivsB.t,derivsB.didmu,output.t);
output.B.dEdmu = spline(derivsB.t,derivsB.dEdmu,output.t);%
output.B.dLdmu = spline(derivsB.t,derivsB.dLdmu,output.t);%
output.B.dQdmu = spline(derivsB.t,derivsB.dQdmu,output.t);%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               plot results                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
set(gca,'fontsize',15,'yscale','log','box','on','NextPlot','add');
warning off;  % this just stops the warning about ignoring negative data on the log scale
semilogy(output.t,abs(output.A.dedmu-output.B.dedmu)./abs(output.A.dedmu),'r--');
semilogy(output.t,abs(output.A.dpdmu-output.B.dpdmu)./abs(output.A.dpdmu),'b--');
semilogy(output.t,abs(output.A.didmu-output.B.didmu)./abs(output.A.didmu),'g--');
semilogy(output.t,abs(output.A.dEdmu-output.B.dEdmu)./abs(output.A.dEdmu),'r');
semilogy(output.t,abs(output.A.dLdmu-output.B.dLdmu)./abs(output.A.dLdmu),'g');
semilogy(output.t,abs(output.A.dQdmu-output.B.dQdmu)./abs(output.A.dQdmu),'b');
warning on; % turn the ability to warn back on
xlabel('t (seconds)');
ylabel('fractional changes in dx/d(ln \mu)');
title(['\epsilon = ' num2str(eps) ',  \delta = ' num2str(delta)]);
legend('x = e','x = p','x = \iota','x = E','x = L','x = Q',4);
set(gca,'ytick',10.^[-14:0]);
axis([-1.5e7   1.2e7   1e-14   1]);



