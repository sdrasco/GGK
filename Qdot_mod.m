function result = Qdot_mod(p,iota,e,Q)
% Eq.60 from Gair&Glampedakis, gr-qc/0510129
%rem Matlab is case sensitive, so q and Q are different
global M spin m;
a = spin;
q = a/M;
q2 = q*q;
s = M/p;
cosi = cos(iota);
sini = sin(iota);

result = sqrt(Q)*(Qdot_over_sqrtQ_2pn(p,iota,e) - ((1-e*e)^1.5)*Qdot_over_sqrtQ_2pn(p,iota,0) ...
    +2.0*((1-e*e)^1.5)*tan(iota)*(Ldot_fit(p,iota,e) + sqrtQ_oversin2i_iotadot_fit(p,iota,e) ));

%disp('stopping here');
%disp('stopped');
