function snr = radarequation(kT,losses,NF,powerTransmit,gainTransmit,gainReceive,receiverBandwidth,lamda)
% RADAREQUATION computes snr for transmitted signal for bistatic SAR
% kT Boltzman's contant times temperature
% losses loss in dB
% noise figure in dB
% transmit power
% gainTransmit transmitter gain
% gainRecieve receiver gain 
snr = (10^((gainTransmit)/10) * 10^((gainReceive)/10)*powerTransmit *lamda^2)./...
    ( (4 * pi)^2 *  10^(NF/10) * kT  * 10^(losses/10) * receiverBandwidth);
