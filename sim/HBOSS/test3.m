%% plot power
% gmtiscript(operating_frequency,SNR,ncpi,Ptransmit,...
%                                                                   radar_height,minmaxrange,width,height,NF,losses,...
%                                                                   T,rcs,azimuth_search,spoil_parameter,...
%                                                                   min_detectable_velocity)
uav = true;
count = 1;
freqs = [5.5 10 17] * 1E9;
if(~uav )
    cpi = [3 5 7];
    snr = [13.5 11.5 10.5];
else
    cpi = [1 3 3];
    snr = [16.5 13.5 13.5];
end

for ii = 1:3
    if(~uav )
        [doppler_dwell_time power_dwell_time] = gmtiscript(freqs(ii),snr(ii),cpi(ii),[],[],[],[],[],[],[],[],[],[],[],[],1);
        
    else
        [doppler_dwell_time power_dwell_time] = gmtiscript(freqs(ii),snr(ii),cpi(ii),150,5,[7 100],.5,.2,[],[],[],[],[],[],[],1);
        
    end
    ppt(count) = sum(power_dwell_time);
    ddt(count) = sum(doppler_dwell_time);
    count = count + 1;


end
count = 1;
freqs = [5.5 10 17] * 1E9;
cpi = [1 1 1];
snr = [16.5 16.5 16.5 ];
for ii = 1:3
    if(~uav )
        [doppler_dwell_time power_dwell_time] = gmtiscript(freqs(ii),snr(ii),cpi(ii),[],[],[],[],[],[],[],[],[],[],[],[],1);
    else
        [doppler_dwell_time power_dwell_time] = gmtiscript(freqs(ii),snr(ii),cpi(ii),150,5,[7 100],.5,.2,[],[],[],[],[],[],[],1);
    end
        sppt(count) = sum(power_dwell_time);
        sddt(count) = sum(doppler_dwell_time);
        count = count + 1;
   

end
%%

att =  max([ppt;ddt]);
satt =  max([sppt;sddt]);
hold off
plot(freqs(1)/1E9,att(1),'x',freqs(2)/1E9,att(2),'o',freqs(3)/1E9,att(3),'d',freqs/1E9,satt,'c+-',freqs/1E9,ppt,'r--',freqs/1E9,ddt,'b--',freqs/1E9,att,'g','markersize',12,'linewidth',1)
if(~uav )
    legend('5.5 GHz, 13.5 dB SNR 15 SIM RECV, 3 CPI','10 GHz, 11.5 dB SNR, 15 SIM RECV, 5 CPI','17 GHz, 10.5 dB SNR, 15 SIM RECV, 7 CPI','STAR LEVEL 2 BOUND 4 SIM RECV','POWER BOUND','DOPPLER BOUND','location','northeast')
else
    legend('5.5 GHz, 16.5 dB SNR, 1 CPI','10 GHz, 13.5 dB SNR, 3 CPI','17 GHz, 13.5 dB SNR, 3 CPI','STAR LEVEL 2 BOUND','POWER BOUND','DOPPLER BOUND','location','northeast')
end
hold on
plot(freqs/1E9,att,'g','linewidth',2)
grid on
xlabel('FREQUENCY GHz')
ylabel('TOTAL DWELL TIME')
if(~uav )
    label = sprintf('C:/figures/time_multi_beam_isr');
else
    label = sprintf('C:/figures/time_multi_beam_uav');
end

print( gcf, '-djpeg', label)