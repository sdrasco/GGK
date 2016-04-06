%Function for comparing my evolution with one example sent by Gair.
global M;
global spin; 
global m;
M = 2.0302544683e5;
spin = 0.9*M;
%just for testing:
%spin = 0.7*M;
%m = 1.0
%
%m = 10.0;
m = 2.0302544683e0;
E0 = m*9.7942675920e-01;
L0 = m*M*2.5407874757e+00
Q0 = m*m*M*M*5.4989046135e+00
%now check that these values correspond to Gair's initial e, p/M, iota
rt_chk = rp_ra(E0,L0,Q0);
ra = rt_chk(1)
rp = rt_chk(2)
p = 2*ra*rp/(ra+rp);
p_over_M = p/M
e = (ra-rp)/(ra+rp)
iota = atan2(sqrt(Q0),abs(L0))

%my unit of tile is solar masses, but Jon's is secs, so have to convert
AU = 499.00478370;
solar_mass = (AU^3.0)*(0.01720209895/86400)^2.0
%above from my unpublished "timing.pdf" 
%t0= 0.e0;
%above E0,L0,Q0 were from 1st row in Jon's tables, which is at 1000s, not
%zero.
t0 = 1.e3/solar_mass
tf = 1.9159e+07/solar_mass

%reading in Jon's soln files
readfileA = fopen('./InspiralELzQ-gair.dat','r');
A = fscanf(readfileA, '%f', [4,inf]);
tg = A(1,:);
Eg = A(2,:);
Lg = A(3,:);
Qg = A(4,:);
%g is for Gair
readfileB = fopen('./Inspiralepiota-gair.dat','r');
B = fscanf(readfileB, '%f', [4,inf]);
eg = B(2,:);
pg = B(3,:);
iotag = acos(B(4,:));
%a mu;
%Nt = [1:Nsteps];
%tspan = t0*ones(1,Nsteps) + ((tf-t0)/(Nsteps -1))*Nt;
%y0 = [E0 L0 Q0];
options = odeset('RelTol', 1e-12);
[trange, solution] = ode45(@ELQdots,[t0 tf],[E0 L0 Q0], options);
Npts = size(solution);
%disp('Npts(1)')
%disp(Npts(1))
%now converting trange to sec to compare with Jon:
trange = trange*solar_mass;
disp('trange');
disp(trange);
disp('solution');
disp(solution);
%just for testing:
%Npts(1) = Npts(1) -3;

%below are quantities that I plot to help debug:
pt = zeros(Npts(1),1);
et = zeros(Npts(1),1);
rat = zeros(Npts(1),1);
rpt = zeros(Npts(1),1);
iotat = zeros(Npts(1),1); 
Edott = zeros(Npts(1),1); 
Ldott = zeros(Npts(1),1);
Qdott = zeros(Npts(1),1);
%
for i =1:Npts(1)
  Et = solution(i,1);
  Lt = solution(i,2); 
  Qt = solution(i,3);
  rt = rp_ra(Et,Lt,Qt);
  ra = rt(1);
  rp = rt(2);
  rat(i) = ra;
  rpt(i) = rp;
  if (ra < rp)
    disp('wrong root order');
    disp('stopped');
  end

  if (ra ~= real(ra)) || (rp ~= real(rp) )
    disp('complex roots');
    disp('stopped');
  end
  %pt(i) = 2*ra*rp/(ra+rp);
  %converting pt to units of M
  pt(i) = 2*ra*rp/(M*(ra+rp));
  et(i) = (ra-rp)/(ra+rp);
  iotat(i) = atan2(sqrt(Qt),abs(Lt));
 % Edott(i) = Edot_mod(pt(i),iotat(i),et(i));
 % Ldott(i) = Ldot_2pn(pt(i),iotat(i),et(i));
 % Qdott(i) = sqrt(solution(i,3))*Qdot_over_sqrtQ_2pn(pt(i),iotat(i),et(i));
  
end
%figure
%plot(solution(:,1)./solution(:,2))
%figure
%plot(trange,Edott)
%figure
%plot(trange,Ldott)
%figure
%plot(trange,Qdott)
%figure
%plot(trange,rat)
%figure
%plot(trange,rpt)
figure
plot(trange,solution(:,1)/m,tg,Eg)
figure
plot(trange,solution(:,2)/(m*M),tg,Lg)
figure
plot(trange,solution(:,3)/(m*m*M*M),tg,Qg)
figure
plot(trange,et,tg,eg)
figure
plot(trange,pt,tg,pg)
figure
plot(trange, iotat,tg,iotag)
%plot(iotat,et)
%plot(pt/M, iotat)
%disp('stopping here');
%disp('stopped');


