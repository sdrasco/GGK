function m = innerprod(gIf,hIf,gIIf,hIIf,delt,Ntot)
%returns the noise-weighted inner product <g,h> = sum of contribution from both channels I and II 
%
%Ntot is number of points in time domain;  in f-space we restrict to
%positive frequencies, so waveforms are half that long.
delf = 1.0/(Ntot*delt);
f = [0:Ntot/2-1]*delf;
SNf = SNlisaCC(f,1.5);
%rem the factor 1.5 above is kappa/T in inverse years

%re-weight g:
gIf = gIf./SNf;
gIIf = gIIf./SNf;
%figure
%plot(cumsum(gIf.*conj(hIf) + hIf.*conj(gIf)));
%figure
%plot(cumsum(gIIf.*conj(hIIf) + hIIf.*conj(gIIf)));
m = 2.0*delf*real(gIf*hIf' + hIf*gIf' + gIIf*hIIf' + hIIf*gIIf');
%disp('stopped');