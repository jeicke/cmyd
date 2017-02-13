function snr=bidetectionrange(Pt,Gt,Gr,lamd,NF,sigma,loss,dwellduty,rng1,rng2);
%snr=bidetectionrange(Txpow (dBw),Txgain (dBi),Rxgain (dBi),wavelength (m),noisefigure (dB),sigma (dBsm),losses (dB),dwell*duty (s),range_1 (m),range_2 (m))
G=Pt+Gt+Gr+20*log10(lamd)-33+228.6-24.6-NF-loss;

snr=G-20*log10(rng1)-20*log10(rng2)+sigma+10*log10(dwellduty);
