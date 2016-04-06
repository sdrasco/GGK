function [DhIf DhIIf] = GenericDerivs(S)
% 
% help is on the way
%


% make sure that the input structure has all the necessary fields
if ~isfield(S,'a') 
    error('GenericDerivs: ''a'' field missing from input structure');
end
if ~isfield(S,'e0') 
    error('GenericDerivs: ''e0'' field missing from input structure');
end
if ~isfield(S,'p0') 
    error('GenericDerivs: ''p0'' field missing from input structure');
end
if ~isfield(S,'r0') 
    error('GenericDerivs: ''r0'' field missing from input structure');
end
if ~isfield(S,'iota0_deg') 
    error('GenericDerivs: ''iota0_deg'' field missing from input structure');
end
if ~isfield(S,'theta0_deg') 
    error('GenericDerivs: ''theta0_deg'' field missing from input structure');
end
if ~isfield(S,'phi0_deg') 
    error('GenericDerivs: ''phi0_deg'' field missing from input structure');
end
if ~isfield(S,'M') 
    error('GenericDerivs: ''M'' field missing from input structure');
end
if ~isfield(S,'mu') 
    error('GenericDerivs: ''mu'' field missing from input structure');
end
if ~isfield(S,'SmallSteps') 
    error('GenericDerivs: ''SmallSteps'' field missing from input structure');
end
if ~isfield(S,'order') 
    error('GenericDerivs: ''order'' field missing from input structure');
end
if ~isfield(S,'thetasb_deg') 
    error('GenericDerivs: ''thetasb_deg'' field missing from input structure');
end
if ~isfield(S,'phisb_deg') 
    error('GenericDerivs: ''phisb_deg'' field missing from input structure');
end
if ~isfield(S,'theta_k_deg') 
    error('GenericDerivs: ''theta_k_deg'' field missing from input structure');
end
if ~isfield(S,'phi_k_deg') 
    error('GenericDerivs: ''phi_k_deg'' field missing from input structure');
end
if ~isfield(S,'D') 
    error('GenericDerivs: ''D'' field missing from input structure');
end
if ~isfield(S,'ParameterKey') 
    error('GenericDerivs: ''ParameterKey'' field missing from input structure');
end
if ~isfield(S,'nparam') 
    error('GenericDerivs: ''nparam'' field missing from input structure');
end
if ~isfield(S,'thetaE_deg') 
    error('GenericDerivs: ''thetaE_deg'' field missing from input structure');
end
if ~isfield(S,'phiE_deg') 
    error('GenericDerivs: ''phiE_deg'' field missing from input structure');
end

% initialize derivative matrices
DhIf = zeros(S.nparam,S.SmallSteps/2); 
DhIIf = zeros(S.nparam,S.SmallSteps/2);

% used to determine finite displacement used in numerical calc of 1st
% derivs
BigEps = 1e-3;
MediumEps = 1e-5;
SmallEps = 1e-8;

% Some of the derivatives are with respect to parameters that are not 
% inputs of the costly InitializeWaveform routine.  We call it with 
% unchanged parameters only once.
W = InitializeWaveform(S);

% initialize extrinsic flag
IsExtrinsic = 0;

% main loop
for i = 1:S.nparam
    
  % send a status report to the screen
  disp(['Calculating derivative with respect to parameter ' ...
      num2str(i) '/' num2str(S.nparam) ': ' S.ParameterKey{i} ' ...']);
  
  % D block
  if strcmp(S.ParameterKey{i},'ln D')
      
      % D-derivative can be done analytically
      [hplus hcross]=ObserveWaveform(W, S.D, S.thetaE_deg, S.phiE_deg, S.order);
      [hI hII tt] = LISAResponse(W, S, hplus, hcross);
      t0 = min(tt);
      delt = abs(tt(1) - tt(2));
      display(['(delt-S.delt) / delt = ' num2str((delt-S.delt) / delt)]); % TEST
      DhIf(i,:) = -cfft(hI,t0,delt,S.SmallSteps);
      DhIIf(i,:) = -cfft(hII,t0,delt,S.SmallSteps);
	  clear hI hII hplus hcross;
      
  % M block    
  elseif strcmp(S.ParameterKey{i},'ln M')
      
      % set extrinsic flag
      IsExtrinsic = 0;
      
      % set eps
      eps = MediumEps; 
      
      % build plus structure
      Sp = S;
      Sp.M = S.M*(1.0+eps);
      
      % build minus structure
      Sm = S;
      Sm.M = S.M*(1.0-eps);
      
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.M - Sm.M)/S.M;
      
  % mu block
  elseif strcmp(S.ParameterKey{i},'ln mu')
      
      % set extrinsic flag
      IsExtrinsic = 0;
      
      % set eps
      eps = SmallEps;

      % build plus structure
      Sp = S;
      Sp.mu = S.mu*(1.0+eps);
      
      % build minus structure
      Sm = S;
      Sm.mu = S.mu*(1.0-eps);
      
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.mu - Sm.mu)/S.mu;

  % a block
  elseif strcmp(S.ParameterKey{i},'a')

      % set extrinsic flag
      IsExtrinsic = 0;
      
      % set eps
      eps = SmallEps;

      % build plus structure
      Sp = S;
      Sp.a = S.a+eps;
      
      % minus
      Sm = S;
      Sm.a = S.a-eps;
      
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.a - Sm.a);

  % p0 block
  elseif strcmp(S.ParameterKey{i},'p0')

      % set extrinsic flag
      IsExtrinsic = 0;
      
      % set eps
      eps = SmallEps;

      % build plus structure
      Sp = S;
      Sp.p0 = S.p0+eps;
      
      % build minus structure
      Sm = S;
      Sm.p0 = S.p0-eps;
      
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.p0 - Sm.p0);
      
  % r0 block    
  elseif strcmp(S.ParameterKey{i},'r0')
      
      % set extrinsic flag
      IsExtrinsic = 0;
      
      % set eps
      eps = SmallEps;

      % build plus structure
      Sp = S;
      Sp.r0 = S.r0+eps;
      
      % build minus structure
      Sm = S;
      Sm.r0 = S.r0-eps;
            
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.r0 - Sm.r0);
      
  % e0 block
  elseif strcmp(S.ParameterKey{i},'e0')
      
      % set extrinsic flag
      IsExtrinsic = 0;
      
      % set eps
      eps = MediumEps; 

      % build plus structure
      Sp = S;
      Sp.e0 = S.e0+eps;
      
      % build minus structure
      Sm = S;
      Sm.e0 = S.e0-eps;
      
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.e0 - Sm.e0);

  % iota0 block
  elseif strcmp(S.ParameterKey{i},'iota0')
      
      % set extrinsic flag
      IsExtrinsic = 0;
      
      % set eps
      eps = SmallEps;

      % build plus structure
      Sp = S;
      Sp.iota0_deg = S.iota0_deg + eps*180/pi;
      
      % build minus structure
      Sm = S;
      Sm.iota0_deg = S.iota0_deg - eps*180/pi;
            
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.iota0_deg - Sm.iota0_deg) * (pi/180);
      
  % theta0 block
  elseif strcmp(S.ParameterKey{i},'theta0')
      
      % set extrinsic flag
      IsExtrinsic = 0;
      
      % set eps
      eps = BigEps;

      % build plus structure
      Sp = S;
      Sp.theta0_deg = S.theta0_deg + eps*180/pi;
            
      % build minus structure
      Sm = S;
      Sm.theta0_deg = S.theta0_deg - eps*180/pi;
            
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.theta0_deg - Sm.theta0_deg) * (pi/180);
      
  % phi0 block    
  elseif strcmp(S.ParameterKey{i},'phi0')
      
      % set extrinsic flag
      IsExtrinsic = 0;
      
      % set eps
      eps = BigEps;

      % build plus structure
      Sp = S;
      Sp.phi0_deg = S.phi0_deg + eps*180/pi;
      
      % build minus structure
      Sm = S;
      Sm.phi0_deg = S.phi0_deg - eps*180/pi;
      
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.phi0_deg - Sm.phi0_deg) * (pi/180);
      
  % thetasb block
  elseif strcmp(S.ParameterKey{i},'thetasb')
      
      % set extrinsic flag
      IsExtrinsic = 1;
      
      % set eps
      eps = BigEps;

      % build plus structure
      Sp = S;
      Sp.thetasb_deg = S.thetasb_deg + eps*180/pi;
      [Sp.thetaE_deg Sp.phiE_deg] = EarthSkyPosition(Sp);
      
      % build minus structure
      Sm = S;
      Sm.thetasb_deg = S.thetasb_deg - eps*180/pi;
      [Sm.thetaE_deg Sm.phiE_deg] = EarthSkyPosition(Sm);
      
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.thetasb_deg - Sm.thetasb_deg) * (pi/180);
      
 % phisb block     
  elseif strcmp(S.ParameterKey{i},'phisb')
      
      % set extrinsic flag
      IsExtrinsic = 1;
      
      % set eps
      eps = BigEps;

      % build plus structure
      Sp = S;
      Sp.phisb_deg = S.phisb_deg + eps*180/pi;
      [Sp.thetaE_deg Sp.phiE_deg] = EarthSkyPosition(Sp);
      
      % build minus structure
      Sm = S;
      Sm.phisb_deg = S.phisb_deg - eps*180/pi;
      [Sm.thetaE_deg Sm.phiE_deg] = EarthSkyPosition(Sm);
      
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.phisb_deg - Sm.phisb_deg) * (pi/180);
      
  % theta_k block
  elseif strcmp(S.ParameterKey{i},'theta_k')
      
      % set extrinsic flag
      IsExtrinsic = 1;
      
      % set eps
      eps = BigEps;

      % build plus structure
      Sp = S;
      Sp.theta_k_deg = S.theta_k_deg + eps*180/pi;
      [Sp.thetaE_deg Sp.phiE_deg] = EarthSkyPosition(Sp);
      
      % build minus structure
      Sm = S;
      Sm.theta_k_deg = S.theta_k_deg - eps*180/pi;
      [Sm.thetaE_deg Sm.phiE_deg] = EarthSkyPosition(Sm);
      
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.theta_k_deg - Sm.theta_k_deg) * (pi/180);
      
  % phi_k block    
  elseif strcmp(S.ParameterKey{i},'phi_k')
      
      % set extrinsic flag
      IsExtrinsic = 1;
      
      % set eps
      eps = BigEps;

      % build plus structure
      Sp = S;
      Sp.phi_k_deg = S.phi_k_deg + eps*180/pi;
      [Sp.thetaE_deg Sp.phiE_deg] = EarthSkyPosition(Sp);
      
      % build minus structure
      Sm = S;
      Sm.phi_k_deg = S.phi_k_deg - eps*180/pi;
      [Sm.thetaE_deg Sm.phiE_deg] = EarthSkyPosition(Sm);
      
      % just in case eps isn't well represented memory
      TrueEps = 0.5*(Sp.phi_k_deg - Sm.phi_k_deg) * (pi/180);
      
  end
  
  % For all parameters other than distance D, the derivatives have to be
  % done numerically.
  if ~strcmp(S.ParameterKey{i},'ln D') 
      
      % plus waveform and response
      [Sp.thetaE_deg Sp.phiE_deg] = EarthSkyPosition(Sp);
      if IsExtrinsic
          [hplusp hcrossp]=ObserveWaveform(W, Sp.D, Sp.thetaE_deg, ...
            Sp.phiE_deg, Sp.order);
          [hIp hIIp tp] = LISAResponse(W, Sp, hplusp, hcrossp);
      else
          Wp = InitializeWaveform(Sp);
          [hplusp hcrossp]=ObserveWaveform(Wp, Sp.D, Sp.thetaE_deg, ...
            Sp.phiE_deg, Sp.order);
          [hIp hIIp tp] = LISAResponse(Wp, Sp, hplusp, hcrossp);
          clear Wp;
      end    
      %clear hplusp hcrossp Sp;
      
      % minus waveform and response
      if IsExtrinsic
          [hplusm hcrossm]=ObserveWaveform(W, Sm.D, Sm.thetaE_deg, ...
            Sm.phiE_deg, Sm.order);
          [hIm hIIm tm] = LISAResponse(W, Sm, hplusm, hcrossm);
      else
          Wm = InitializeWaveform(Sm);
          [hplusm hcrossm]=ObserveWaveform(Wm, Sm.D, Sm.thetaE_deg, ...
            Sm.phiE_deg, Sm.order);
          [hIm hIIm tm] = LISAResponse(Wm, Sm, hplusm, hcrossm);
          clear Wm;          
      end
      %clear hplusm hcrossm Sm;
      
      % reset extrinsic flag
      IsExtrinsic = 0;
      
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
 
      % derivatives of hI and hII
      dhI = (hIp - hIm)/(2.0*TrueEps);
      %dhI = (hIp - hIm)/(2.0*eps);
      clear hIp hIm;
      dhII = (hIIp - hIIm)/(2.0*TrueEps);
      %dhII = (hIIp - hIIm)/(2.0*eps);
      clear hIIp hIIm;
      
      % Fourier-transform the derivatives
      t0 = min(tt);
      delt = abs(tt(1) - tt(2));
      display(['(delt-S.delt) / delt = ' num2str((delt-S.delt) / delt)]); % TEST
      DhIf(i,:) = cfft(dhI,t0,delt,S.SmallSteps);
      DhIIf(i,:) = cfft(dhII,t0,delt,S.SmallSteps);       
      clear dhI dhII tt;
      
  end
  
end


end
