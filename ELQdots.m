function dy = ELQdots(t,y)
% ELQdots  the rhs of the numerical kludge ODEs for E(t),L_z(t),Q(t)
%       based on the modified 2PN formulae in Gair Glampedakis, gr-qc/0510129
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
if (ra < rp)
    disp('wrong root order');
    disp('stopped');
end

if (ra ~= real(ra)) || (rp ~= real(rp) )
    disp('complex roots');
    disp('stopped');
end
p = 2*ra*rp/(ra+rp);
e = (ra-rp)/(ra+rp);
iota = atan2(sqrt(Q),abs(L));
q = a/M;
q2 = q*q;
s = M/p;
cosi = cos(iota);
sini = sin(iota);

dy(1) = Edot_mod(p,iota,e);   %eq.20, with GHK replaced by "2PN" for Edot terms and "mod" for Ldot and Qdot

dy(2) = Ldot_mod(p,iota,e);    %eqs. 59,45, and 57

dy(3) = Qdot_mod(p,iota,e,Q);%eqs. 60, 56, 57 and 58

%dy(1) = Edot_2pn(p,iota,e);   %eq.20, with GHK replaced by "2PN" for Edot terms and "mod" for Ldot and Qdot

%dy(1) = Edot_mod(p,iota,e); 

%dy(2) = Ldot_2pn(p,iota,e);    %eqs. 59,45, and 57

%dy(3) = sqrt(Q)* Qdot_over_sqrtQ_2pn(p,iota,e);

dy = dy';




