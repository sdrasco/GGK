function r = rp_ra(E,L,Q)
%used to find perihelion and aphelion for any (E,L,Q), by finding roots of Eq. 1.
%Has been tested using Drasco code that makes the reverse transformation.
global M spin m;
a = spin;
%from eq. 6 of gair et al
c4 = E^2 - m^2;
c3 = 2.0*M*m^2;
c2 = 2.0*E^2*a^2 - 2.0*a*L*E - (L - a*E)^2 -Q - m^2*a^2;
c1 = 2.0*M*((L - a*E)^2 + Q);
%c0 = E^2*a^4 - 2*a^3*L*E + a^2*L^2 - a^2*((L - a*E)^2 + Q);
%note almost all terms in c0 cancel, except the last one.  all that
%remains is:
c0 = -Q*a^2;
p = [c4 c3 c2 c1 c0];

% find roots analytically
r = roots(p);

% zero any complex roots
r(find(~isreal(r))) = 0;

% note: Often the first root r(1) is ra (the larger one), and r(2) 
% is rp.  However, roots() doesn't return things sorted by size, so reorder
% them by hand.  This is important near merger.
r = flip(sort(r));

% if we don't have two real roots, we faild.
if r(1) == 0 || r(2) == 0
   error('rp_ra(): unable to find real values of rp and ra.');
end


end

