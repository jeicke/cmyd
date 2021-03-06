function r=detectionrange(Pt,Gt,Gr,lamd,NF,sigma,loss,dwellduty,SNR);
%r=detectionrange(Txpow,Txgain,Rxgain,wavelength,noisefigure,sigma,losses,dwell*duty,requiredsnr)
G=Pt+Gt+Gr+20*log10(lamd)-33+228.6-24.6-NF-loss;

r4=G-SNR+sigma+10*log10(dwellduty);
r=10.^(r4/40);