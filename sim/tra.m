parameters_real;
compute_stap = true;
velocity = 270;
pulses = 64;
radar_velocity = [velocity*sind(45);velocity * cosd(45); 0];
Npulses = pulses ;
make_targets = false;
make_clutter = true;
make_noise = false;
transmitters = 1;
resolution = .125/2;
prf = 3000;
mimo;
imagesc(f_simo,ranges/1000,20*log10(abs(fftshift( Dbf_simo,2))))
title('SIMO NO RANGE AMG');

caxis([0 70])
colorbar
drawnow;
xlabel('FREQUENCY (HZ)')
ylabel('RANGE (KM)')
print('-djpeg70','C:/Users/rnl/Desktop/figures/SIMO NO RANGE AMG')
figure;

parameters_real;
range_ambiguous = false;
compute_stap = true;
velocity = 270;
Npulses = 64;
make_targets = false;
make_clutter = true;
make_noise = false;
transmitters = 2;
resolution = .125/2;
prf = 3000;
radar_velocity = [velocity*sind(45);velocity * cosd(45); 0];
mimo;
figure
sinr_stap_simo_u = sinr_stap_simo;
sinr_stap_mimo_u = sinr_stap_mimo;
plot(f_simo,10 * log10(sinr_stap_simo),'b',f_mimo,10 * log10(sinr_stap_mimo)-3,'r','linewidth',2);
grid
xlabel('FREQUENCY (HZ)')
ylabel('SINR (dB)')
legend('SIMO','MIMO');
print('-djpeg70','C:/Users/rnl/Desktop/figures/SINRNRA')
figure;
imagesc(f_mimo,ranges/1000,20*log10(abs(fftshift( Dbf_mimo,2))))
title('MIMO NO RANGE AMG');
caxis([0 70])
colorbar
drawnow;
xlabel('FREQUENCY (HZ)')
ylabel('RANGE (KM)')
print('-djpeg70','C:/Users/rnl/Desktop/figures/MIMO NO RANGE AMG')

figure;
parameters_real;

compute_stap = true;
range_ambiguous = true;
velocity = 270;
Npulses = 64;
make_targets = false;
make_clutter = true;
make_noise = false;
transmitters = 1;
resolution = .125/2;
prf = 3000;
radar_velocity = [velocity*sind(45);velocity * cosd(45); 0];
mimo;
imagesc(f_simo,unambiguous_range/1000,20*log10(abs(fftshift( Dbf_simo,2))))
title('SIMO RANGE AMG');
caxis([0 70])
colorbar
drawnow;
xlabel('FREQUENCY (HZ)')
ylabel('RANGE (KM)')
print('-djpeg70','C:/Users/rnl/Desktop/figures/SIMO RANGE AMG')
figure;

parameters_real;

compute_stap = true;
range_ambiguous = true;
velocity = 270;
Npulses = 64;
make_targets = false;
make_clutter = true;
make_noise = false;
transmitters = 2;
resolution = .125/2;
prf = 3000;
radar_velocity = [velocity*sind(45);velocity * cosd(45); 0];
mimo;
imagesc(f_mimo,unambiguous_range/1000,20*log10(abs(fftshift( Dbf_mimo,2))))
title('MIMO RANGE AMG');
caxis([0 70])
colorbar
drawnow;
xlabel('FREQUENCY (HZ)')
ylabel('RANGE (KM)')
print('-djpeg70','C:/Users/rnl/Desktop/figures/MIMO RANGE AMG')
figure
plot(f_simo,10 * log10(sinr_stap_simo),'b',f_mimo,10 * log10(sinr_stap_mimo)-3,'r',...
     f_simo,10 * log10(sinr_stap_simo_u),'b--',f_mimo,10 * log10(sinr_stap_mimo_u)-3,'r--','linewidth',2);
grid
xlabel('FREQUENCY (HZ)')
ylabel('SINR (dB)')
legend('SIMO','MIMO','SIMO NRA','MIMO NRA');
print('-djpeg70','C:/Users/rnl/Desktop/figures/SINRALL')
