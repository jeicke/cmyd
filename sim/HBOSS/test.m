%% plot power
% gmtiscript(operating_frequency,SNR,ncpi,Ptransmit,...
%                                                                   radar_height,minmaxrange,width,height,NF,losses,...
%                                                                   T,rcs,azimuth_search,spoil_parameter,...
%                                                                   min_detectable_velocity)
uav = false;
count = 1;
freqs = [5.5 8 10 15 17] * 1E9;
snrs = [10.5 13.5 16.5];
for sn = snrs
    sn
    count2 = 1;
    for tf = freqs
        if(uav)
            [doppler_dwell_time power_dwell_time] = gmtiscript(tf,sn,[],150,5,[7 100],.6,.6,[],[],[],[],[],[],[],1);
        else
            [doppler_dwell_time power_dwell_time] = gmtiscript(tf,sn,[],[],[],[],[],[],[],[],[],[],[],[],[],1);
        end
        ppt(count,count2) = sum(power_dwell_time);
        ddt(count,count2) = sum(doppler_dwell_time);
        count2 = count2 + 1;
    end
    count = count + 1;
end
%%
plot(freqs/1E9,ppt(1,:),'r-x',freqs/1E9,ppt(2,:),'b-o',freqs/1E9,ppt(3,:),'g-d','linewidth',2)
xlabel('FREQUENCY GHz')
ylabel('TOTAL DWELL TIME FOR SNR')
legend('10.5 dB SNR','13.5 dB SNR','16.5 dB SNR')
grid on
if(uav)
    label = sprintf('C:/figures/power_summary_uav');
else
    label = sprintf('C:/figures/power_summary_isr');
end
print( gcf, '-djpeg', label )
%figure
%plot(freqs/1E9,ddt(1,:),'r-x',freqs/1E9,ddt(2,:),'b-o',freqs/1E9,ddt(3,:),'g-d','linewidth',2)
%xlabel('FREQUENCY GHz')
%ylabel('TOTAL DWELL TIME FOR 1 M/S OF DOPPLER RESOLUTION')
%grid on
%label = sprintf('C:/figures/power_summary_isr');
%print( gcf, '-djpeg', label )