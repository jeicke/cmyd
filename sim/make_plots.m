%% plot of beampattern on target for mimo and simo
parameters;
az = [-10:.125:10];
Npulses = 64;
make_targets = true;
make_clutter = false;
make_noise = false;
transmitters = 4;
resolution = .125;
prf = 1500;
mimo;

Npulses = 64;
prf = 1500;
transmitters = 1;
mimo;
p_mimo = 20 * log10(abs(squeeze(Dbf_sweep_mimo(342,1,:))));
p_simo = 20 * log10(abs(squeeze(Dbf_sweep_simo(342,1,:))));
plot(az,p_mimo,'b',az,p_simo,'r','linewidth',2)

legend('MIMO','SIMO')
grid on;
xlabel('ANGLE (THETA)')
ylabel('POWER TARGET(dB)')
%print('-djpeg70','C:/Users/rnl/Desktop/figures/azimuthal_beampattern')
return
%% plot of elevation beampattern on target for mimo and simo
figure;
parameters;

el = [-20:.5:20];
Npulses = 128;
make_targets = true;
make_clutter = false;
make_noise = false;
transmitters = 1;
resolution = .125;
prf = 1500;
mimo;
Npulses = 128;
prf = 1500;
transmitters = 2;
mimo;
p_mimo = 20 * log10(abs(squeeze(Dbf_sweep_el_mimo(342,1,:))));
p_simo = 20 * log10(abs(squeeze(Dbf_sweep_el_simo(342,1,:))));
plot(el,p_mimo,'b',el,p_simo,'r','linewidth',2)

legend('MIMO','SIMO')
grid on;
xlabel('ANGLE (ELEVATION) (phi)')
ylabel('POWER TARGET(dB)')
print('-djpeg70','C:/Users/rnl/Desktop/figures/elevation_beampattern')
%% azimuthal clutter beampattern
figure;
clear Dbf_sweep_mimo;
clear Dbf_sweep_simo;
parameters;
az = [-10:.125:10];
Npulses = 64;
make_targets = false;
make_clutter = true;
make_noise = true;
transmitters = 1;
resolution = .125/2;
prf = 1500;
mimo;
Npulses = 64;
prf = 1500;
transmitters = 2;
mimo;
p_mimo = 20 * log10(sum(abs(squeeze(Dbf_sweep_mimo(320:360,1,:))),1));
p_simo = 20 * log10(sum(abs(squeeze(Dbf_sweep_simo(320:360,1,:))),1));
plot(az,p_mimo,'b',az,p_simo,'r','linewidth',2)
legend('MIMO','SIMO')
grid on;
xlabel('ANGLE (THETA)')
ylabel('POWER CLUTTER(dB)')
print('-djpeg70','C:/Users/rnl/Desktop/figures/azimuthal_beampattern_clutter')
%% width of clutter null based on theoretical beam patterns
% HGA.theta_source = 0;
% HGA.phi_source = EL(320) * 180/pi;
% HGA.SLx = mimo_receive_geometry;
% HGA.subarray_geometry = subarray_geometry;
% HGA.subarray_window = subarray_window ;
% [transmit_beam_pattern_mimo] = calcBeamPattern(operating_frequency,AZ * 180/pi,EL(320) * 180/pi,[],HGA);
% HGA.SLx = fullarray_geometry;
% [transmit_beam_pattern_simo] = calcBeamPattern(operating_frequency,AZ * 180/pi,EL(320) * 180/pi,[],HGA);
% f_mimo
% dop_mimo = interp1(dopplers(320,:),transmit_beam_pattern,f_mimo,'linear',0);
%% various plots
v = 18;
f = 5E9;
c = 299792458;
lamda = c/f;
w = 1.2;
EL = 0;
theta = 0;
vangle = 0:.1:90;
format  bank
df1 =  2 * v/lamda * cosd(EL) * 2 * (cosd(vangle-theta)-cosd(vangle + asind(lamda/w) - theta));
w = 2 * w;
df2 =  2 * v/lamda * cosd(EL) * 2 * (cosd(vangle-theta)-cosd(vangle + asind(lamda/w) - theta));
w12 = df1-df2;

w = 2.4;
EL = 0;
theta = 0;
vangle = 0:.1:90;
format  bank
df1 =  2 * v/lamda * cosd(EL) * 2 * (cosd(vangle-theta)-cosd(vangle + asind(lamda/w) - theta));
w = 2 * w;
df2 =  2 * v/lamda * cosd(EL) * 2 * (cosd(vangle-theta)-cosd(vangle + asind(lamda/w) - theta));
w24 = df1-df2;
plot(vangle,w12,'r',vangle,w24,'g','linewidth',2)
grid
xlabel('ANGLE BETWEEN LOOK AND MOTION')
ylabel('(CLUTTER WIDTH SIMO - CLUTTER WIDTH MIMO) Hz')
legend('WIDTH 1.2 m V 18 m/s','WIDTH 2.4m V 18 m/s','location','northwest');
print('-djpeg70','C:/Users/rnl/Desktop/figures/clutterwidthdifferencevslookangle')





w = 2.4;
EL = -90:0;
theta = 0;
vangle = 90;
format  bank
df1 =  2 * v/lamda * cosd(EL) * 2 * (cosd(vangle-theta)-cosd(vangle + asind(lamda/w) - theta));
w = 2 * w;
df2 =  2 * v/lamda * cosd(EL) * 2 * (cosd(vangle-theta)-cosd(vangle + asind(lamda/w) - theta));
w24 = df1-df2;
figure;
plot(EL,w24,'g','linewidth',2)
grid
xlabel('ELEVATION ANGLE')
ylabel('(CLUTTER WIDTH SIMO - CLUTTER WIDTH MIMO) Hz')
legend('WIDTH 2.4 m V 18 m/s','location','northwest');
print('-djpeg70','C:/Users/rnl/Desktop/figures/clutterwidthdifferencevselevationangle')


%% width of clutter null for various dangles
figure;
parameters;
velocity = 80;
pulses = 256;
radar_velocity = [velocity; 0; 0];
Npulses = pulses ;
make_targets = false;
make_clutter = true;
make_noise = true;
transmitters = 1;
resolution = .125/4;
prf = 1500;
mimo;
Npulses = pulses ;
prf = 1500;
transmitters = 2;
mimo;
plot_beam(f_simo,f_mimo,Dbf_mimo,Dbf_simo,sum(abs(D_mimo),3),sum(abs(D_simo(:,:,:)),3),320:360)
title('90 degrees')
print('-djpeg70','C:/Users/rnl/Desktop/figures/90degreeclutterwidth')
figure;
parameters;
radar_velocity = [0; velocity; 0];
Npulses = pulses ;
make_targets = false;
make_clutter = true;
make_noise = true;
transmitters = 1;
resolution = .125/4;
prf = 1500;
mimo;
Npulses = pulses ;
prf = 1500;
transmitters = 2;
mimo;
plot_beam(f_simo,f_mimo,Dbf_mimo,Dbf_simo,sum(abs(D_mimo),3),sum(abs(D_simo(:,:,:)),3),320:360)
title('0 degrees')
print('-djpeg70','C:/Users/rnl/Desktop/figures/0degreeclutterwidth')
figure;
parameters;
radar_velocity = [velocity*sind(45);velocity * cosd(45); 0];
Npulses = pulses ;
make_targets = false;
make_clutter = true;
make_noise = true;
transmitters = 1;
resolution = .125/4;
prf = 1500;
mimo;
Npulses = pulses ;
prf = 1500;
transmitters = 2;
mimo;
plot_beam(f_simo,f_mimo,Dbf_mimo,Dbf_simo,sum(abs(D_mimo),3),sum(abs(D_simo(:,:,:)),3),320:360)
title('45 degrees')
print('-djpeg70','C:/Users/rnl/Desktop/figures/45degreeclutterwidth')
%% calculate SINR for pre-doppler STAP
figure;
parameters;
compute_stap = true;
velocity = 120;
pulses = 128;
radar_velocity = [velocity*sind(45);velocity * cosd(45); 0];
Npulses = pulses ;
make_targets = false;
make_clutter = true;
make_noise = false;
transmitters = 1;
resolution = .125/2;
prf = 1500;
mimo;
Npulses = pulses ;
prf = 1500;
transmitters = 2;
mimo;
sinr_stap_simo = mean(sinr_stap_simo,1);
sinr_stap_mimo = mean(sinr_stap_mimo,1);
plot(f_simo,10 * log10(sinr_stap_simo),'b',f_mimo,10 * log10(sinr_stap_mimo)-3,'r','linewidth',2);
grid
xlabel('FREQUENCY (HZ)')
ylabel('SINR (dB)')
legend('SIMO','MIMO');
print('-djpeg70','C:/Users/rnl/Desktop/figures/sinrplot')
%% calculate SINR for post-doppler STAP
figure;
parameters;
compute_stap = true;
velocity = 18;
pulses = 256;
radar_velocity = [velocity*sind(45);velocity * cosd(45); 0];
Npulses = pulses ;
make_targets = false;
make_clutter = true;
make_noise = false;
transmitters = 1;
resolution = .125/2;
prf = 1500;
mimo;
Npulses = pulses ;
prf = 1500;
transmitters = 2;
mimo;
s1 = fftshift(mean(sinr_stap_simo(ii-10:ii+10,:),1));
s2 = fftshift(mean(sinr_stap_mimo(ii-10:ii+10,:),1));
plot(f_simo,10 * log10(s1),'b',f_mimo,10 * log10(s2),'r','linewidth',2);
%axis([f_simo(1) f_simo(end) -100 2])
grid
xlabel('FREQUENCY (HZ)')
ylabel('SINR (dB)')
legend('SIMO','MIMO');
print('-djpeg70','C:/Users/rnl/Desktop/figures/sinrplot')