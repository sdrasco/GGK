function result = sqrtQ_oversin2i_iotadot_fit(p,iota,e)
% Eq.58 from Gair&Glampedakis, gr-qc/0510129, multiplied by sqrt{Q}/sin^2i
%
global M spin m;
a = spin;
q = a/M;
q2 = q*q;
q3= q2*q;
q4 = q3*q;
s = M/p;
sqrt_s = s^0.5;
cosi = cos(iota);
sini = sin(iota);
cos2i = cosi*cosi;
cos3i = cos2i*cosi;
cos4i = cos3i*cosi;
cos5i = cos4i*cosi;
d1a = -10.742;  d1b = 28.5942; d1c = -9.07738; d2a = -1.42836; d2b = 10.7003;
d2c = -33.709; c1a = -28.1517; c1b = 60.9607; c1c = 40.9998; c2a = -0.348161;
c2b = 2.37258; c2c = -66.6584; c3a = -0.715392; c3b = 3.21593; c3c = 5.28888;
c4a = -7.61034; c4b = 128.878; c4c = -475.465; c5a = 12.2908; c5b = -113.125;
c5c = 306.119; c6a = 40.9259; c6b = -347.271; c6c = 886.503; c7a = -25.4831;
c7b = 224.227; c7c = -490.982; c8a = -9.00634; c8b = 91.1767; c8c = -297.002;
c9a = -0.645; c9b = -5.13592; c9c = 47.1982; f1a = -283.955; f1b = 736.209;
f2a = 483.266; f2b = -1325.19; f3a = -219.224; f3b = 634.499; f4a = -25.8203;
f4b = 82.078; f5a = 301.478; f5b = -904.161; f6a = -271.966; f6b = 827.319;

c10a = -0.0309341; c10b = -22.2416; c10c = 7.55265; c11a = - 3.33476; c11b = 22.7013;
c11c = -12.47; f7a = -162.268; f7b = 247.168; f8a = 152.125; f8b = -182.165;
f9a = 184.465; f9b = -267.553; f10a = -188.132; f10b = 254.067;


result = 6.4*(m*m/M)*q*(s^5.0)*( (61.e0/24.e0) + s*(d1a + d1b*sqrt_s + d1c*s) +q2*s*(d2a + d2b*sqrt_s + d2c*s) ...
    + q*cosi*sqrt_s*(c10a + c10b*s + c10c*s^1.5) + q2*cos2i*s*(c11a + c11b*sqrt_s + c11c*s) ...
    + q3*cosi*(s^2.5)*(f7a + f7b*sqrt_s + q*(f8a +f8b*sqrt_s) + cos2i*(f9a + f9b*sqrt_s) +q*cos2i*(f10a + f10b*sqrt_s)));

%disp('stopping here');
%disp('stopped');
