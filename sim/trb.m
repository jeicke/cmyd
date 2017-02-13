figure;
parameters3;
compute_stap = true;
velocity = 18;
pulses = 64;
radar_velocity = [velocity*sind(5);velocity * cosd(5); 0];
Npulses = pulses ;
make_targets = false;
make_clutter = true;
make_noise = false;
transmitters = 1;
resolution = .125;
prf = 1500;
mimo;
Npulses = pulses ;
prf = 1500;
transmitters = 4;
mimo;
sinr_stap_simoU = mean(sinr_stap_simo,1);
sinr_stap_mimoU = mean(sinr_stap_mimo,1);
plot(f_simo,10 * log10(sinr_stap_simoU),'b',f_mimo,10 * log10(sinr_stap_mimoU)-3,'r','linewidth',2);
grid
xlabel('FREQUENCY (HZ)')
ylabel('SINR (dB)')
legend('SIMO','MIMO');
print('-djpeg70','C:/Users/rnl/Desktop/figures/sinrplot')

figure;
parameters3;
range_ambiguous = true;
compute_stap = true;
velocity = 18;
pulses = 64;
radar_velocity = [velocity*sind(5);velocity * cosd(5); 0];
Npulses = pulses ;
make_targets = false;
make_clutter = true;
make_noise = false;
transmitters = 1;
resolution = .125;
prf = 1500;
mimo;
Npulses = pulses ;
prf = 1500;
transmitters = 4;
mimo;
plot(f_simo,10 * log10(sinr_stap_simo),'b',f_mimo,10 * log10(sinr_stap_mimo)-3,'r',...
    'linewidth',2);
grid
xlabel('FREQUENCY (HZ)')
ylabel('SINR (dB)')
legend('SIMO','MIMO');
print('-djpeg70','C:/Users/rnl/Desktop/figures/sinrplotboth')

figure;
sinr_stap_simo = mean(sinr_stap_simo,1);
sinr_stap_mimo = mean(sinr_stap_mimo,1);
sinr_stap_mimoU = mean(sinr_stap_mimoU,1);
sinr_stap_simoU = mean(sinr_stap_simoU,1);
plot(f_simo,10 * log10(sinr_stap_simo),'b',f_mimo,10 * log10(sinr_stap_mimo),'r',...
    f_simo,10 * log10(sinr_stap_simoU),'b--',f_mimo,10 * log10(sinr_stap_mimoU),'r--');
grid
xlabel('FREQUENCY (HZ)')
ylabel('SINR (dB)')
legend('SIMO','MIMO','SIMO NO RANGE AMBIGUOUS','MIMO NO RANGE AMBIGUOUS');
print('-djpeg70','C:/Users/rnl/Desktop/figures/sinrplotra')