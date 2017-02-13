function Pt=radarpower(Gt,Ar,SNR,NF,sigma,loss,dwellduty,rng)
%Pt=radarpower(Tx Antenna Gain (dB), Rx Antenna Area (m^2),SNR (dB),NF dB),sigma (dB),loss (dB),dwellduty ,rng (m))
G=Gt+10*log10(Ar)-22+228.6-24.6-NF-loss;
Pt=SNR+40*log10(rng)-G-sigma-10*log10(dwellduty);
