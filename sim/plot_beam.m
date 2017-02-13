function plot_beam(f_simo,f_mimo,bmimo,bsimo,ubmimo,ubsimo,range_gates)

p_mimo = fftshift(20 * log10(mean(abs(bmimo(range_gates,:)),1)));
p_simo = fftshift(20 * log10(mean(abs(bsimo(range_gates,:)),1)));

p_ubmimo = fftshift(20 * log10(mean(abs(ubmimo(range_gates,:)),1)));
p_ubsimo = fftshift(20 * log10(mean(abs(ubsimo(range_gates,:)),1)));
plot(f_simo,p_simo,'b',f_simo,p_ubsimo,'b--',f_mimo,p_mimo,'r',f_mimo,p_ubmimo,'r--','linewidth',2)
grid
legend('SIMO','ubsimo','MIMO','ubmimo','location','west')
xlabel('FREQUENCY (Hz)')
ylabel('POWER (dB)')