function dy = ELQdots(t,y)
% ELQdots the rhs of the numerical kludge ODEs for E(t),L_z(t),Q(t)
% based on the formulae in GG = Gair Glampedakis, gr-qc/0510129
%
% These are modified 2PN flux rules (modified to better handle small
% eccentricty, large inclination, and to better fit Hughe's circular
% Teukolsky code).
%
global M spin m;
a = spin;
dy = [0,0,0];
E = y(1);
L = y(2); 
Q = y(3);
r = rp_ra(E,L,Q);
ra = r(1);
rp = r(2);

% debugging: show me rp, ra, and p as time evolves
%plot(t,ra/M,'ko',t,rp/M,'kx',t,(2*ra*rp/(ra+rp))/M,'r+');
%drawnow

p = 2*ra*rp/(ra+rp); % note: p has dimensions
e = (ra-rp)/(ra+rp);
iota = atan2(sqrt(Q),abs(L));
q = a/M;
q2 = q*q;
s = M/p;
cosi = cos(iota);
sini = sin(iota);

% this block gives best agreement with what plots in GG
% Edot: GG's eq.20, with Edot_GHK replaced by eq 44 and "mod" for Ldot and Qdot
dy(1) = Edot_mod(p,iota,e);  

% Ldot: GG's eqs. 59,45, and 57
dy(2) = Ldot_mod(p,iota,e);  

% Qdot: GG's eqs. 60, 56, 57 and 58
dy(3) = Qdot_mod(p,iota,e,Q);

dy = dy';




