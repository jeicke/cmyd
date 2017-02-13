%% test mimo RA
figure;
parameters;
radar_position= [0; 0; 10];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf = 1500;
Npulses =50;
range_ambiguous = false;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 550;
compute_radarequation =true;
mimo;
sinr_simo = 10*log10(sinr);

parameters;
v = 100;
radar_position= [0; 0; 10];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =1500;
range_ambiguous = false;
transmitters = 2;
angle = 60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =50;
range_gate = 550;
compute_radarequation =true;
mimo;
sinr_mimo = 10*log10(sinr);

parameters;
radar_position= [0; 0; 10];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =1500;
Npulses =50;
range_ambiguous = true;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 550;
compute_radarequation =true;
mimo;
sinr_simo_ra = 10*log10(sinr);

parameters;
v = 100;
radar_position= [0; 0; 10];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =1500;
angle =  60;
range_ambiguous = true;
transmitters = 2;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =50;
range_gate = 550;
compute_radarequation =true;
mimo;
sinr_mimo_ra = 10*log10(sinr);

transmitters = 4;
Npulses = 30;
%mimo;
f = doppler_bins';
sinr_mimo2 = nan * ones(size(f));%10*log10(sinr);

plot(f,sinr_simo,'r',...
    f,sinr_simo_ra,'g',...
    f,(sinr_mimo),'k',...
    f,(sinr_mimo_ra),'b',...
    f,sinr_mimo2-6,'g','linewidth',2)
grid on;
legend('SIMO','SIMO RANGE AMBIGUOUS','MIMO','MIMO RANGE AMBIGUOUS','location','southwest')
xlabel('DOPPLER')
ylabel('SINR (dB)');
axis([f(1) f(end) -59 5])
param_text = ['MIMO RANGE AMBIGUOUS'  sprintf('\nPRF %d HZ PULSES %d VELOCITY %d',prf,Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/MIMO_RANGE_AMBIGUOUS0')