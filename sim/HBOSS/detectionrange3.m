function snr=detectionrange(Pt,Gt,Gr,lamd,NF,sigma,loss,dwellduty,rng);
%r=detectionrange(Txpow,Txgain,Rxgain,wavelength,noisefigure,sigma,losses,dwell*duty,requiredsnr)
G=Pt+Gt+Gr+20*log10(lamd)-33+228.6-24.6-NF-loss;
snr=G-40*log10(rng)+sigma+10*log10(dwellduty);
