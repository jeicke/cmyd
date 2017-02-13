hold off;
plot(y,10*log10(abs(sar)))
grid on
 xlabel('CROSS RANGE (m)')
 ylabel('POWER (dB)')
yy = [40 10 0 -20 -50 -80];
mmin = 40;
mmax = 100;
hold on;
t = M/prf
dy = t * 100
bw = (10000*lamda/(dy))
for ii = 1:length(yy)
    plot([yy(ii) yy(ii)],[mmin mmax],'r','linewidth',.5)
   % plot([yy(ii)-bw/2 yy(ii)-bw/2],[mmin mmax],'r--','linewidth',.5)
   % plot([yy(ii)+bw/2 yy(ii)+bw/2],[mmin mmax],'r--','linewidth',.5)
end
%axis([-100 100 60 90])