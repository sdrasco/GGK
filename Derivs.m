function [DhIf DhIIf] = Derivs(nparam, a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
        sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k)

%1: ln D
%2: ln M
%3: ln mu
%4: a \equiv S/M^2
%5: p0 \equiv p/M at time t0
%6  r0 \equiv r/M at time t0
%7: e0
%8: iota0 in radians
%9: theta_0 in radians
%10: phi0 (radians) -- this is azimuthal location of particle (wrt to MBH spin) at t=t0,
%where phi =0 is azimuthal direction from source to Earth
%11: thetasb in radians
%12:phisb (radians)
%13:theta_k in radians
%14:phi_k (radians)

DhIf = zeros(14,SmallSteps/2); 
DhIIf = zeros(14,SmallSteps/2);
Cdelt = (tF -tI)/(SmallSteps -1);

nparam=14;
%eps = 1.e-5;  %This is finite displacement used in numerical calc of 1st derivs; 
%note Curt has NOT yet tried to vary eps and show robustness --needs to be checked.

BigEps = 1e-3;
SmallEps = 1e-8;
%BigEps = 1e-5;
%SmallEps = 1e-5;

for i = 1:nparam
  disp(['i =', num2str(i)]);
  if (i==1)
        [hI hII tp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        delt = abs(tp(1) - tp(2));
        display(['(delt-Cdelt) / Cdelt = ' num2str((delt-Cdelt) / Cdelt)]); % TEST
        tp0 = min(tp);
        DhIf(1,:) = -cfft(hI,tp0,delt,SmallSteps);
        DhIIf(1,:) = -cfft(hII,tp0,delt,SmallSteps);
        clear hI; clear hII;
  end

  if (i > 1)
        
     if (i==2)
        eps = SmallEps;
        Mp = M*(1.0+eps);
        Mm = M*(1.0-eps);
        [hIp hIIp tp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  Mp, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm tm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  Mm, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
    
     elseif (i==3)
        eps = SmallEps;
        mup = mu*(1.0+eps);
        mum = mu*(1.0-eps);
        [hIp hIIp tp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mup, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm tm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mum, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        
     elseif (i==4)
        eps = SmallEps;
        ap = a+eps;
        am = a-eps;
        [hIp hIIp tp] = LISAh2(ap, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm tm] = LISAh2(am, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
      
     elseif (i==5)
        eps = SmallEps;
        p0p = p0 + eps;
        p0m = p0 - eps;
        [hIp hIIp tp] = LISAh2(a, e0, p0p, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm tm] = LISAh2(a, e0, p0m, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        
     elseif (i==6)
        eps = SmallEps;
        r0p = r0 + eps;
        r0m = r0 - eps;
        [hIp hIIp tp] = LISAh2(a, e0, p0, iota_deg0, r0p, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm tm] = LISAh2(a, e0, p0, iota_deg0, r0m, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
      
     elseif (i==7)
        eps = SmallEps;
        e0p = e0 + eps;
        e0m = e0 - eps;
        [hIp hIIp tp] = LISAh2(a, e0p, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm tm] = LISAh2(a, e0m, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
     
     elseif (i==8)
        eps = SmallEps;
        iota_deg0p = iota_deg0 + eps*180.0/pi;
        iota_deg0m = iota_deg0 - eps*180.0/pi;
        [hIp hIIp tp] = LISAh2(a, e0, p0, iota_deg0p, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm tm] = LISAh2(a, e0, p0, iota_deg0m, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);

     elseif (i==9)
        eps = BigEps;
        theta0_degp = theta0_deg + eps*180.0/pi;
        theta0_degm = theta0_deg - eps*180.0/pi;
        [hIp hIIp tp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_degp, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm tm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_degm, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);

     elseif (i==10)
        eps = BigEps;
        phi0_degp = phi0_deg + eps*180.0/pi;
        phi0_degm = phi0_deg - eps*180.0/pi;
        [hIp hIIp tp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_degp, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm tm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_degm, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
		  
     elseif (i==11)
        eps = BigEps;
        thetasbp = thetasb + eps;
        thetasbm = thetasb - eps;
        [hIp hIIp tp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasbp, phisb, theta_k, phi_k);
        [hIm hIIm tm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasbm, phisb, theta_k, phi_k);
		  
	elseif (i==12)
        eps = BigEps;
        phisbp = phisb + eps;
        phisbm = phisb - eps;
        [hIp hIIp tp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisbp, theta_k, phi_k);
        [hIm hIIm tm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisbm, theta_k, phi_k);

    elseif (i==13)
        eps = BigEps;
        theta_kp = theta_k + eps;
        theta_km = theta_k - eps;
        [hIp hIIp tp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_kp, phi_k);
        [hIm hIIm tm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_km, phi_k);
		  
    elseif (i==14)
        eps = BigEps;
        phi_kp = phi_k + eps;
        phi_km = phi_k - eps;
        [hIp hIIp tp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_kp);
        [hIm hIIm tm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_km);
     end
     
     % before taking differences of hI and hII, we interpolate the plus
     % and minus waveform functions in order to get them at equivalent
     % times
     tmin = max([min(tp) min(tm)]);
     tmax = min([max(tp) max(tm)]);
     tt = linspace(tmin,tmax,length(tp));
     hIp = spline(tp,hIp,tt);
     hIIp = spline(tp,hIIp,tt);
     hIm = spline(tm,hIm,tt);
     hIIm = spline(tm,hIIm,tt);
     clear tp tm;

     % if eps is too small, something has gone wrong
     if eps < 1e-14
         error(['Derivs: eps has been set to a bad value: eps = ' num2str(eps)]);
     end
     
     % get timestep and initial time for FFT
     delt = abs(tt(1) - tt(2));
     tt0 = min(tt);
     display(['(delt-Cdelt) / Cdelt = ' num2str((delt-Cdelt) / Cdelt)]); % TEST
     
     dhI = (hIp - hIm)/(2.0*eps);
     dhII = (hIIp - hIIm)/(2.0*eps);
     DhIf(i,:) = cfft(dhI,tt0,delt,SmallSteps);
     DhIIf(i,:) = cfft(dhII,tt0,delt,SmallSteps);
     
     % triggers an error if I forgot to set eps somewhere
     clear eps;
     
     %now making plot to check derivs
     %tt = 0:SmallSteps-1;
     %tt = tI + tt*(tF-tI)/(SmallSteps-1);
     %SmallStepsI = floor(SmallSteps*(t0-tI)/(tF-tI)); 
     %fig_null = makeplots(tt,SmallStepsI,dhI,dhII);
       
     clear hIp; clear hIIp; clear hIm; clear hIIm; clear dhI; clear dhII
  end
end

%disp('stopping here for debugging purposes');
%disp('stopped');
