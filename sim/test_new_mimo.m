%% test of SINR as a function of angle
parameters;
Npulses = 64;
v = 100;
angle =  0:30:90;
count = 1;
clear sinr_all;
for ii = 0:30:90
    radar_velocity = v * [sind(ii); cosd(ii); 0];
    mimo;
    sinr_all(:,count) = sinr;
    count = count + 1;
end
f = doppler_bins';%prf/transmitters * (-Npulses/2:Npulses/2-1)/(Npulses/2);

sinr_all = 10*log10(sinr_all);
plot(f,sinr_all(:,1),'r',...
    f,sinr_all(:,2),'c',...
    f,sinr_all(:,3),'b',...
    f,sinr_all(:,4),'k','linewidth',3)

hold on;
% plot theory curves
color = {'r--','c--','b--','k--'}
color2 = {'r-.','c-.','b-.','k-.  '}
count = 1;
for ii = 0:30:90
    s = mod(2 * v/(c/operating_frequency) * cosd(ii) * cos(EL(range_gate)),prf);
    if(s>prf/2)
        s = s - prf;
    end
    df =s
    plot([df df],[-30 5],color{count},'linewidth',2)
    
    w = width;
    
    df1 =  2 * v/(c/operating_frequency) * cos(EL(range_gate)) * 2 * (cosd(ii)-cosd(ii +  17/20 * 1/2 * asind((c/operating_frequency)/w)))
    plot([df-df1/2 df-df1/2],[-30 5], color2{count},[df+df1/2 df+df1/2],[-30 5], color2{count},'linewidth',2);
     count = count + 1;
end


grid on;
legend('0','30','60','90','location','southwest')
xlabel('DOPPLER (Hz)')
ylabel('SINR (dB)');
axis([f(1) f(end) -30 5])
param_text = ['SINR VERSUS ANGLE ' sprintf('\nPRF %d HZ PULSES %d VELOCITY %d',prf,Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/SINR_VERSUS_ANGLE')
%% test range folding
figure;
parameters;
v= 100;
count = 1;
clear sinr_all;
clear sinr_all_ra

for ii = [0 60 90]
    parameters;
     
    radar_velocity =v * [sind(ii); cosd(ii); 0];
    mimo;
    sinr_all(:,count) = sinr;
    range_ambiguous = true;
   
    mimo;
    sinr_all_ra(:,count) = sinr;
    count = count + 1;
end
f = doppler_bins';
sinr_all = 10*log10(sinr_all);
sinr_all_ra = 10*log10(sinr_all_ra);
plot(f,sinr_all(:,1),'r',...
    f,sinr_all_ra(:,1),'r--',...
    f,sinr_all(:,2),'c',...
    f,sinr_all_ra(:,2),'c--',...
    f,sinr_all(:,3),'k',...
    f,sinr_all_ra(:,3),'k--','linewidth',2)
grid on;
legend('0','0 RANGE AMBIGUOUS','60','60 RANGE AMBIGUOUS','90','90 RANGE AMBIGUOUS','location','southwest')
title(param_text);
xlabel('DOPPLER')
ylabel('SINR (dB)');
axis([f(1) f(end) -30 5])
param_text = ['SINR SHOWING EFFECTS OF RANGE AMBIGUOUS CLUTTER ' sprintf('\nPRF %d HZ PULSES %d VELOCITY %d',prf,Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/SINR_WITH_RAC')
%% test range folding vs PRF
figure;
parameters;
v= 100;
count = 1;
clear sinr_all;
clear sinr_all_ra

clear f;
clear sinr_all_ra;
clear sinr_all;
for ii = [300 1000 1500]
    parameters;
   
    range_gate = 70;    
    prf = ii;

    radar_velocity =v * [sind(60); cosd(60); 0];
    mimo;
    sinr_all{:,count} = sinr;
    range_ambiguous = true;
   
    mimo;
    sinr_all_ra{:,count} = sinr;
    f{:,count} = doppler_bins';
    count = count + 1;
end

plot(f{:,1},10*log10(sinr_all{:,1}),'r',...
    f{:,1},10*log10(sinr_all_ra{:,1}),'r--',...
    f{:,2},10*log10(sinr_all{:,2}),'c',...
    f{:,2},10*log10(sinr_all_ra{:,2}),'c--',...
    f{:,3},10*log10(sinr_all{:,3}),'k',...
    f{:,3},10*log10(sinr_all_ra{:,3}),'k--','linewidth',2)
grid on;
legend('300','300 RANGE AMBIGUOUS','1000','1000 RANGE AMBIGUOUS','1500','1500 RANGE AMBIGUOUS','location','southwest')
xlabel('DOPPLER')
ylabel('SINR (dB)');
axis([-750 750 -30 5])
param_text = ['RANGE FOLDING VS PRF' sprintf('\nPULSES %d VELOCITY %d',Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/SINR_VS_PRF')
%% test range folding vs PRF with range attenuation
figure;
parameters;
v= 100;
count = 1;
clear sinr_all;
clear sinr_all_ra
ii = 60;
clear f;
clear sinr_all_ra;
clear sinr_all;

for ii = [1000 1500]
    parameters;
    Npulses = 50;
    compute_radarequation =true;
    transmit_beam = [0 -4.9];
    receive_beam  = [0 -4.9];
    range_gate = 340;    
    prf = ii;

    radar_velocity =v * [sind(60); cosd(60); 0];
    mimo;
    sinr_all{:,count} = sinr;
    range_ambiguous = true;
   
    mimo;
    sinr_all_ra{:,count} = sinr;
    f{:,count} = doppler_bins';
    count = count + 1;
end

plot(f{:,1},10*log10(sinr_all{:,1}),'r',...
    f{:,1},10*log10(sinr_all_ra{:,1}),'r--',...
    f{:,2},10*log10(sinr_all{:,2}),'c',...
    f{:,2},10*log10(sinr_all_ra{:,2}),'c--','linewidth',2)
grid on;
legend('1000','1000 RANGE AMBIGUOUS','1500','1500 RANGE AMBIGUOUS','location','southwest')
xlabel('DOPPLER')
ylabel('SINR (dB)');
axis([-750 750 -60 5])
param_text = ['RANGE FOLDING WITH RADAR EQUATION' sprintf('\nPULSES %d VELOCITY %d',Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/RANGEFOLDING_RADAREQUATION')
%% test range attenuation

figure;
parameters;
v = 100;
count = 1;
clear sinr_all;
clear sinr_all_ra;
clear sinr_all_att;
f = doppler_bins';
for ii = [60]
    parameters;
    power_clipping = 40;
    range_gate = 120;
    radar_velocity =v* [sind(ii); cosd(ii); 0];
     range_ambiguous = true;
    compute_radarequation =true;
    mimo;
    sinr_all_att(:,count) = sinr;
    
     range_ambiguous = false;
    compute_radarequation = false;
    
    mimo;
    sinr_all(:,count) = sinr;
    range_ambiguous = true;
    
    mimo;
    sinr_all_ra(:,count) = sinr;
    count = count + 1;
end
sinr_all = 10*log10(sinr_all);
sinr_all_ra = 10*log10(sinr_all_ra);
sinr_all_att = 10 * log10(sinr_all_att);
plot(f,sinr_all(:,1),'r',...
    f,sinr_all_ra(:,1),'k',...
    f,sinr_all_att(:,1),'b',...
    'linewidth',2)
grid on;
legend('60','60 RANGE AMBIGUOUS','60 RANGE ATTENUATED','location','southwest')
xlabel('DOPPLER')
axis([f(1) f(end) -30 5])
param_text = ['EFFECT OF RADAR EQUATION'  sprintf('\nPRF %d HZ PULSES %d VELOCITY %d',prf,Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/RADAREQUATION')

%% test 2-d array
figure;
parameters;
v = 100;
count = 1;
clear sinr_all;
clear sinr_all_ra;
clear sinr_all_att;
f = doppler_bins';
for ii = [60]
     
    parameters;
    power_clipping = 40;
    radar_velocity = v * [sind(ii); cosd(ii); 0];
    range_ambiguous = false;
    compute_radarequation = false;
    line_array = false;
    compute_subarray = true;
    range_gate = 200;
    transmit_beam = [0 -9.92];
    receive_beam  = [0 -9.92];
    mimo;
    sinr_all(:,count) = sinr;
    f = doppler_bins;
    
    parameters;
    power_clipping = 40;
    radar_velocity = v * [sind(ii); cosd(ii); 0];
    range_ambiguous = true;
    compute_radarequation = true;
    compute_subarray = true;
    line_array =false;
    range_gate = 200;
    transmit_beam = [0 -9.92];
    receive_beam  = [0 -9.92];
    mimo;
    sinr_all_att(:,count) = sinr;
    f_att = doppler_bins;
    
    range_ambiguous = true;
    line_array = false;
    compute_subarray = true;
    range_gate = 200;
    transmit_beam = [0 -9.92];
    receive_beam  = [0 -9.92];
    mimo;
    sinr_all_ra(:,count) = sinr;
    f_ra = doppler_bins;
    count = count + 1;
end
sinr_all = 10*log10(sinr_all);
sinr_all_ra = 10*log10(sinr_all_ra);
sinr_all_att = 10 * log10(sinr_all_att);
plot(f,sinr_all(:,1),'r',...
    f_ra,sinr_all_ra(:,1),'k',...
    f_att,sinr_all_att(:,1),'b',...
    'linewidth',2)
grid on;
legend('60','60 RANGE AMBIGUOUS','60 RANGE ATTENUATED','location','southwest')
xlabel('DOPPLER')
ylabel('SINR (dB)');
axis([f(1) f(end) -59 5])
param_text = ['2-D ARRAY'  sprintf('\nPRF %d HZ PULSES %d VELOCITY %d',prf,Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/2DARRAY')

%% test mimo
close all;
clear all;
parameters;
Npulses =59;
range_ambiguous = false;
v = 200;
angle =  90;
radar_velocity = v * [sind(angle); cosd(angle); 0];
width = .6;
mimo;
sinr_simo = 10*log10(sinr);
f = doppler_bins';

transmitters = 2;
Npulses =60;
width = .6;
mimo;
sinr_mimo = 10*log10(sinr);
f1 = doppler_bins';

transmitters = 4;
Npulses = 60;
width = .6;
mimo;
sinr_mimo2 = 10*log10(sinr);
f2 = doppler_bins';

plot(f,sinr_simo,'r',...
    f1,sinr_mimo,'k',...
    f2,sinr_mimo2,'g','linewidth',2)
grid on;
legend('SIMO','MIMO','MIMO*2','location','southwest')
xlabel('DOPPLER')
ylabel('SINR (dB)');
param_text = ['MIMO'  sprintf('\nPRF %d HZ PULSES %d VELOCITY %d',prf,Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/MIMO')

%axis([f(1) f(end) -59 5])
%% test mimo RA with radar equation
parameters;
radar_position= [0; 0; .25];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf = 1500;
Npulses =50;
range_ambiguous = false;
v = 100;
angle =  60;
compute_radarequation =true;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 60;
mimo;
sinr_simo = 10*log10(sinr);

parameters;
v = 100;
radar_position= [0; 0; .25];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =1500;
range_ambiguous = false;
transmitters = 2;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =50;
range_gate = 60;
compute_radarequation =true;
mimo;
sinr_mimo = 10*log10(sinr);

parameters;
radar_position= [0; 0; .25];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =1500;
Npulses =50;
range_ambiguous = true;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 60;
compute_radarequation =true;
mimo;
sinr_simo_ra = 10*log10(sinr);

parameters;
v = 100;
radar_position= [0; 0; .25];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =1500;
angle =  60;
range_ambiguous = true;
transmitters = 2;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =50;
range_gate = 60;
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
    f,sinr_mimo2,'g','linewidth',2)
grid on;
legend('SIMO','SIMO RANGE AMBIGUOUS','MIMO','MIMO RANGE AMBIGUOUS','location','southwest')
xlabel('DOPPLER')
ylabel('SINR (dB)');
axis([f(1) f(end) -59 5])
param_text = ['MIMO RANGE AMBIGUOUS WITH RADAR EQUATION'  sprintf('\nPRF %d HZ PULSES %d VELOCITY %d',prf,Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/MIMO_RANGE_AMBIGUOUS_RADAR_EQUATION')
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
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/MIMO_RANGE_AMBIGUOUS_ALT10')
%% test mimo RA 2
parameters;
line_array = true;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf = 400;
Npulses =50;
range_ambiguous = false;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 340;
compute_radarequation =false;
mimo;
sinr_simo = 10*log10(sinr);

parameters;
v = 100;
line_array = true;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
range_ambiguous = false;
transmitters = 2;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =50;
range_gate = 340;
mimo;
sinr_mimo = 10*log10(sinr);

parameters;
line_array = true;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
Npulses =50;
range_ambiguous = true;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 340;
mimo;
sinr_simo_ra = 10*log10(sinr);

parameters;
line_array = true;
v = 100;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
angle =  60;
range_ambiguous = true;
transmitters = 2;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =50;
range_gate = 340;
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
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/MIMO_RANGE_AMBIGUOUS2')
%% test mimo RA 2 with Range Attenuation
parameters;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf = 400;
Npulses =50;
range_ambiguous = false;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 340;
compute_radarequation =true;
mimo;
sinr_simo = 10*log10(sinr);

parameters;
v = 100;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
range_ambiguous = false;
transmitters = 2;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =50;
range_gate = 340;
compute_radarequation =true;
mimo;
sinr_mimo = 10*log10(sinr);

parameters;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
Npulses =50;
range_ambiguous = true;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 340;
compute_radarequation =true;
mimo;
sinr_simo_ra = 10*log10(sinr);

parameters;
v = 100;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
angle =  60;
range_ambiguous = true;
transmitters = 2;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =50;
range_gate = 340;
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


param_text = ['MIMO RANGE AMBIGUOUS WITH RADAR EQUATION'  sprintf('\nPRF %d HZ PULSES %d VELOCITY %d',prf,Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/MIMO_RANGE_AMBIGUOUS_RADAREQUATION')

axis([f(1) f(end) -59 5])
%% test mimo RA 2 with Range Attenuation close
parameters;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf = 400;
Npulses =50;
range_ambiguous = false;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 140;
compute_radarequation =true;
mimo;
sinr_simo = 10*log10(sinr);

parameters;
v = 100;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
range_ambiguous = false;
transmitters = 2;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =50;
range_gate = 140;
compute_radarequation =true;
mimo;
sinr_mimo = 10*log10(sinr);

parameters;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
Npulses =50;
range_ambiguous = true;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 140;
compute_radarequation =true;
mimo;
sinr_simo_ra = 10*log10(sinr);

parameters;
v = 100;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
angle =  60;
range_ambiguous = true;
transmitters = 2;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =50;
range_gate = 140;
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


param_text = ['MIMO RANGE AMBIGUOUS WITH RADAR EQUATION'  sprintf('\nPRF %d HZ PULSES %d VELOCITY %d',prf,Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/MIMO_RANGE_AMBIGUOUS_RADAREQUATION_CLOSE')

axis([f(1) f(end) -59 5])
%% test pre-doppler stap
parameters;
Npulses =30;
v = 100;
radar_position= [0; 0; 10];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =1500;
range_gate = 550;
range_ambiguous = true;
angle = 60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
transmitters = 1;
compute_stap = false;
compute_radarequation =false;
mimo;
sinr_reg = sinr;

parameters;
Npulses =30;
v = 100;
radar_position= [0; 0; 10];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =1500;
range_gate = 550;
transmitters = 1;
angle = 60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
compute_stap = true;
compute_radarequation =false;
range_ambiguous =true;
mimo;
sinr_pds = sinr;

parameters;
transmitters = 2;
Npulses =30;
v = 100;
radar_position= [0; 0; 10];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =1500;
range_gate = 550;
angle = 60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
compute_stap = true;
compute_radarequation =false;
range_ambiguous = true;


mimo;
sinr_pds_mimo = sinr;

f = doppler_bins';
plot(f,10*log10(sinr_pds),'r',f,10*log10(sinr_reg),'g',f,10*log10(sinr_pds_mimo),'b','linewidth',2)
grid on;
legend('PRE-DOPPLER STAP','STAP','PRE-DOPPLER MIMO','location','southwest')
xlabel('DOPPLER')
ylabel('SINR (dB)');
%axis([f(1) f(end) -59 8])
param_text = ['PRE DOPPLER STAP'  sprintf('\nPRF %d HZ PULSES %d',prf,Npulses)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/PRE_DOPPLER_STAP')
%% test pre-doppler stap with ra 
parameters;

Npulses =3;
v = 100;
angle =  90;
radar_velocity = v * [sind(angle); cosd(angle); 0];

compute_stap = false;
compute_radarequation =false;
mimo;
sinr_reg = sinr;

parameters;

Npulses = 3;
v = 100;
angle =  90;
radar_velocity = v * [sind(angle); cosd(angle); 0];
compute_stap = true;
compute_radarequation =false;
mimo;
sinr_pds = sinr;

parameters;

transmitters = 2;
Npulses = 3;
v = 100;
angle =  90;
radar_velocity = v * [sind(angle); cosd(angle); 0];
compute_stap = true;
compute_radarequation =false;
mimo;
sinr_pds_mimo = sinr;

f = doppler_bins';
plot(f,10*log10(sinr_pds),'r',f,10*log10(sinr_reg),'g',f,10*log10(sinr_pds_mimo),'b','linewidth',2)
grid on;
legend('PRE-DOPPLER STAP','STAP','PRE-DOPPLER MIMO','location','southwest')
xlabel('DOPPLER')
ylabel('SINR (dB)');
axis([f(1) f(end) -59 8])
param_text = ['PRE DOPPLER STAP'  sprintf('\nPRF %d HZ PULSES %d',prf,Npulses)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/PRE_DOPPLER_STAP')
%% test mimo RA 2 with Range Attenuation and pre-doppler stap
parameters;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf = 400;
Npulses =60;
range_ambiguous = false;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 340;
compute_radarequation =true;
compute_stap = true;
mimo;
sinr_simo = 5*log10(sinr);

parameters;
v = 100;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
range_ambiguous = false;
transmitters = 2;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =60;
range_gate = 340;
compute_radarequation =true;
compute_stap = true;
mimo;
sinr_mimo = 5*log10(sinr);

parameters;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
Npulses =60;
range_ambiguous = true;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
range_gate = 340;
compute_radarequation =true;
compute_stap = true;
mimo;
sinr_simo_ra = 5*log10(sinr);

parameters;
v = 100;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf =400;
angle =  60;
range_ambiguous = true;
transmitters = 2;
radar_velocity = v * [sind(angle); cosd(angle); 0];
Npulses =60;
range_gate = 340;
compute_radarequation =true;
compute_stap = true;
mimo;
sinr_mimo_ra = 5*log10(sinr);

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


param_text = ['MIMO RANGE AMBIGUOUS WITH RADAR EQUATION PRE-DOPLER STAP'  sprintf('\nPRF %d HZ PULSES %d VELOCITY %d',prf,Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/MIMO_RANGE_AMBIGUOUS_RADAREQUATION_PREDOPPLER')

axis([f(1) f(end) -59 5])
%% estimate covariance
parameters;
radar_position= [0; 0; 5];
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];
prf = 400;
Npulses =50;
range_ambiguous = false;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];

compute_radarequation =true;
range_gate = [];
mimo;
sinr_simo = 10*log10(sinr);
%% compute post doppler stap
parameters;

v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
compute_post_doppler_stap = true;
mimo;
sinr_postdoppler = 10*log10(sinr);

parameters;
v = 100;
angle =  60;
radar_velocity = v * [sind(angle); cosd(angle); 0];
compute_post_doppler_stap = false;
mimo;
sinr = 10*log10(sinr);
f = doppler_bins';
plot(f,sinr_postdoppler,'r',...
    f,sinr,'g',...
    'linewidth',2)
grid on;
legend('SIMO POST DOPPLER','SIMO','location','southwest')
xlabel('DOPPLER')
ylabel('SINR (dB)');


param_text = ['SIMO POST-DOPLER STAP'  sprintf('\nPRF %d HZ PULSES %d VELOCITY %d',prf,Npulses,v)];
title(param_text);
print('-djpeg70','C:/Documents and Settings/jeicke/My Documents/mimofigures/SIMO_POST_DOPPLER')

axis([f(1) f(end) -59 5])
