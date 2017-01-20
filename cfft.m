function hf = cfft(ht,tI,delt,Ntot)
%uses Matlab fft, but throws away redundant second half (the negative
%frequencies), and uses normalization closer to the continuous
%(non-discrete) fft.

hf = fft(ht);
hf = hf(1:Ntot/2);
delf = 1.0/(Ntot*delt);
f = [0:Ntot/2-1]*delf;
time_shift_fac = exp(i*2.0*pi*f*tI);
hf = delt*hf.*time_shift_fac;
%disp('stopped');
