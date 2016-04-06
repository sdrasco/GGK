function result = Ldot_fit(p,iota,e)
% Eq.57 from Gair&Glampedakis, gr-qc/0510129
%Note Ldot_fit doesn't actually have any dependence on the variable e; it's
%based on fitting to circular orbits in Kerr
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


result = -6.4*(m*m/M)*(s^3.5)*( cosi + q*(s^1.5)*(61.e0/24.e0  - (61.e0/8.e0)*cos2i) ...
    -(1247.e0/336.e0)*s*cosi + 4.e0*pi*(s^1.5)*cosi - (44711.e0/9072.e0)*(s^2.0)*cosi + ...
    q2*(s^2.0)*cosi*(33.e0/16.e0 - (45.0/8.0)*sini*sini ) + (s^2.5)*  ...
    (q*(d1a + d1b*sqrt_s + d1c*s) +q3*(d2a + d2b*sqrt_s + d2c*s) + cosi*(c1a + c1b*sqrt_s +c1c*s) ...
     + q2*cosi*(c2a + c2b*sqrt_s + c2c*s) + q4*cosi*(c3a + c3b*sqrt_s + c3c*s) ...
     +q*cos2i*(c4a + c4b*sqrt_s + c4c*s) + q3*cos2i*(c5a + c5b*sqrt_s + c5c*s) ...
     +q2*cos3i*(c6a + c6b*sqrt_s + c6c*s) + q4*cos3i*(c7a + c7b*sqrt_s + c7c*s) ...
     +q3*cos4i*(c8a + c8b*sqrt_s + c8c*s) + q4*cos5i*(c9a + c9b*sqrt_s + c9c*s) ) ...
     +(s^3.5)*q*cosi*( f1a + f1b*sqrt_s + q*(f2a + f2b*sqrt_s) + q2*(f3a + f3b*sqrt_s) ...
     +cos2i*(f4a + f4b*sqrt_s) + q*cosi*cosi*(f5a + f5b*sqrt_s) + q2*cos2i*(f6a + f6b*sqrt_s)));
     

%disp('stopping here');
%disp('stopped');
