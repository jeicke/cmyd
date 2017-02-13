function dwellduty=detectionrange(Pt,Gt,Gr,lamd,NF,sigma,loss,rng,SNR);
%dwellduty=detectionrange(Txpow,Txgain,Rxgain,wavelength,noisefigure,sigma,losses,range(m),requiredsnr)
G=Pt+Gt+Gr+20*log10(lamd)-33+228.6-24.6-NF-loss;

%r4=G-SNR+sigma+10*log10(dwellduty);
%r=10.^(r4/40);

dwellduty=10*log10(rng.^4)+SNR-sigma-G;
dwellduty=10.^(dwellduty/10);