function [DhIf DhIIf] = CEDerivs(eps,nparam, a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
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

nparam=10;

DhIf = zeros(nparam,SmallSteps/2); 
DhIIf = zeros(nparam,SmallSteps/2);
delt = (tF -tI)/(SmallSteps -1);


%eps = 1.e-5;  %This is finite displacement used in numerical calc of 1st derivs; 
%note Curt has NOT yet tried to vary eps and show robustness --needs to be checked.
%eps = 1.e-7; 

for i = 1:nparam
  disp(['Calculating derivative with respect to parameter ' num2str(i) '/' num2str(nparam) ' ...']);
  if (i==1)
        [hI hII] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        DhIf(1,:) = -cfft(hI,tI,delt,SmallSteps);
        DhIIf(1,:) = -cfft(hII,tI,delt,SmallSteps);
        clear hI; clear hII;
  end

  if (i > 1)
       
     if (i==2)
     elseif (i==2)
        Mp = M*(1.0+eps);
        Mm = M*(1.0-eps);
        TrueEps = 0.5*(Mp - Mm)/M;
        [hIp hIIp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  Mp, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  Mm, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
    
     elseif (i==3)
        mup = mu*(1.0+eps);
        mum = mu*(1.0-eps);
        TrueEps = 0.5*(mup - mum)/mu;
        [hIp hIIp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mup, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mum, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        
     elseif (i==4)
        ap = a+eps;
        am = a-eps;
        TrueEps = 0.5*(ap - am);
        [hIp hIIp] = LISAh2(ap, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm] = LISAh2(am, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
      
     elseif (i==5)
        p0p = p0 + eps;r0p = p0*(1 +e0/3.0); 
        p0m = p0 - eps;r0m = p0*(1 +e0/3.0); 
        TrueEps = 0.5*(p0p - p0m);
        [hIp hIIp] = LISAh2(a, e0, p0p, iota_deg0, r0p, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm] = LISAh2(a, e0, p0m, iota_deg0, r0m, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        
     elseif (i==6)
        phi0_degp = phi0_deg + eps;%*180.0/pi;
        phi0_degm = phi0_deg - eps;%*180.0/pi;
        TrueEps = 0.5*(phi0_degp - phi0_degm);
        [hIp hIIp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_degp, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
        [hIm hIIm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_degm, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_k);
		  
     elseif (i==7)
        thetasbp = thetasb + eps;
        thetasbm = thetasb - eps;
        TrueEps = 0.5*(thetasbp - thetasbm);
        [hIp hIIp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasbp, phisb, theta_k, phi_k);
        [hIm hIIm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasbm, phisb, theta_k, phi_k);
		  
	elseif (i==8)
        phisbp = phisb + eps;
        phisbm = phisb - eps;
        TrueEps = 0.5*(phisbp - phisbm);
        [hIp hIIp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisbp, theta_k, phi_k);
        [hIm hIIm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisbm, theta_k, phi_k);

    elseif (i==9)
        theta_kp = theta_k + eps;
        theta_km = theta_k - eps;
        TrueEps = 0.5*(theta_kp - theta_km);
        [hIp hIIp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_kp, phi_k);
        [hIm hIIm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_km, phi_k);
		  
    elseif (i==10)
        phi_kp = phi_k + eps;
        phi_km = phi_k - eps;
        TrueEps = 0.5*(phi_kp - phi_km);
        [hIp hIIp] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_kp);
        [hIm hIIm] = LISAh2(a, e0, p0, iota_deg0, r0, theta0_deg, phi0_deg, ... 
          sign_rdot0, sign_Tdot0,  M, mu, tI, t0,  tF, SmallSteps, BigSteps, tol, coords, D, thetasb, phisb, theta_k, phi_km);
    end
  
     dhI = (hIp - hIm)/(2*eps);
     dhII = (hIIp - hIIm)/(2*eps);
     %dhI = (hIp - hIm)/(2*TrueEps);
     %dhII = (hIIp - hIIm)/(2*TrueEps);
     DhIf(i,:) = cfft(dhI,tI,delt,SmallSteps);
     DhIIf(i,:) = cfft(dhII,tI,delt,SmallSteps);
     
     %now making plot to check derivs
     %tt = 0:SmallSteps-1;
     %tt = tI + tt*(tF-tI)/(SmallSteps-1);
     %SmallStepsI = floor(SmallSteps*(t0-tI)/(tF-tI)); 
     %fig_null = makeplots(tt,SmallStepsI,dhI,dhII);
       
     clear hIp; clear hIIp; clear hIm; clear hIIm; clear dhI; clear dhII;
  end
end

disp(['Finished calculating derivatives.']);

%disp('stopping here for debugging purposes');
%disp('stopped');
