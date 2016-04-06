% Compares h waveforms calculated by both the time technique and freq technique.

%%%%% Initial Parameters %%%%%

S.a = 0.5;
S.e0 = 0.2;
S.p0 = 15;
S.iota0_deg = 33;
%S.r0 = S.p0 / ( 1 - S.e0);
S.r0 = S.p0 / ( 1 + S.e0);
%S.theta0_deg = 90;
S.theta0_deg = 90 - S.iota0_deg;
S.phi0_deg = 0;
S.sign_rdot0 = 1;
%S.sign_Tdot0 = -1;
S.sign_Tdot0 = 1;
S.M = 1.4*8;
S.mu = 1.4;
S.t0 = 0;
S.tspan = 4000*(S.M*4.9255e-6);
S.BigSteps = 10;
S.SmallSteps = 7000;
S.tol = 1e-7;
S.coords = 'spherical';

% these last parameters tell where where the binary is in the sky, and how
% far it is (S.D)
S.D = 1.e17;
S.thetasb_deg = 45;
S.phisb_deg = 150;
S.theta_k_deg = 60;
S.phi_k_deg = 200;
S.order = 'quadrupole';



%%%%% Time Domain Calc %%%%%

%Make waveform
waveform_time = InitializeWaveform(S);

%Find theta and phi location
[thetaE_deg phiE_deg] = EarthSkyPosition(S);

%Observe waveform
[hplus hcross] = ObserveWaveform(waveform_time, S.D, thetaE_deg, phiE_deg, S.order);%%

h_time = hplus - (1i)*hcross;



% %%%%% Frequency Domain Calc %%%%%
% 
% % initialize the waveform
% waveform_frq = InitializeWaveformSeries(S);
% 
% % spline it
% output = SplineWaveformSeries(waveform_frq);
% 
% %Make time vectors
% tt_in_M = linspace(output.t(1),output.t(end),S.SmallSteps);
% tt_in_sec = tt_in_M * (waveform_frq.M*4.9255e-6);
% 
% %Create structure to input into EvolvingSeriesEval
% hstruc.N = length(output.mkn);
% hstruc.mkn = output.mkn;
% hstruc.coefs_pp = output.coefs_pp;
% hstruc.freq1_pp = output.OmegaPhi_pp;
% hstruc.freq2_pp = output.OmegaTheta_pp;
% hstruc.freq3_pp = output.OmegaR_pp;
% hstruc.t = output.t;
% 
% %Solve series for waveform
% h_frq = EvolvingSeriesEval(hstruc,tt_in_M);



%%%%% Plot Waveforms %%%%%

%figure; plot(tt_in_M,h_time,'b',tt_in_M,h_frq,'r')

