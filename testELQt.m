%Function for testing evolution eqs for p,e,iota. It plots curve
%that SHOULD look like lower solid curve in Fig 8 of Gair&Glamp..,
%gr/qc/051019
global M;
global spin; 
global m;
M = 1.e6;
spin = 0.8*M;
%just for testing:
%spin = 0.7*M;
%m = 1.0
%
m = 10.0;
r = 3.3e1*M;
E0 = m*(1.0 - 0.5*M/r);
L0 = 0.837*m*sqrt(M*r);
Q0 = 1.e-16*L0^2;
%Q0 = 1.e-6*L0*L0;
t0= 0.e0;
tf = 4.5286984272e8*M;
%above values of a,r,L0 give e=0.6 at p=20, a=0.8, and take
%evolution almost to end, so can
%be compared to Fig.8
%Above is correct tf when I use Ldot_mod and Qdot_mod in the Edot_mod
%formula
%
%To see that Ldot_2pn and Qdot_2pn are the wrong things to stick in that
%formula, try it, and try:
%tf = 4.43e8*M
%You'll see (e,p) track curves the wrong way at end.

%a mu;
%Nt = [1:Nsteps];
%tspan = t0*ones(1,Nsteps) + ((tf-t0)/(Nsteps -1))*Nt;
%y0 = [E0 L0 Q0];
options = odeset('RelTol', 1e-8);
[trange, solution] = ode45(@ELQdots,[t0 tf],[E0 L0 Q0], options);
Npts = size(solution);
%disp('Npts(1)')
%disp(Npts(1))
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
  pt(i) = 2*ra*rp/(ra+rp);
  et(i) = (ra-rp)/(ra+rp);
  iotat(i) = atan2(sqrt(Qt),abs(Lt));
  Edott(i) = Edot_mod(pt(i),iotat(i),et(i));
  Ldott(i) = Ldot_2pn(pt(i),iotat(i),et(i));
  Qdott(i) = sqrt(solution(i,3))*Qdot_over_sqrtQ_2pn(pt(i),iotat(i),et(i));
  
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
plot(trange,solution(:,1))
%figure
%plot(trange,solution(:,2))
%figure
%plot(trange,solution(:,3))
%figure
%plot(trange,et)
%figure
%plot(trange,pt)
figure
plot(pt/M,et)
%plot(iotat,et)
%plot(pt/M, iotat)
%disp('stopping here');
%disp('stopped');


