%% plot power
% gmtiscript(operating_frequency,SNR,ncpi,Ptransmit,...
%                                                                   radar_height,minmaxrange,width,height,NF,losses,...
%                                                                   T,rcs,azimuth_search,spoil_parameter,...
%                                                                   min_detectable_velocity)

uav = false;
count = 1;
freqs = [5.5 10 17] * 1E9;
if(~uav )
    cpi = [3 5 7];
    snr = [13.5 11.5 10.5];
    elev_ranges{1} = [300 30 25];
    elev_ranges{2} = [300 48 25];
    elev_ranges{3} = [300 80 48 35 25];
    az_spoil = [1 2 4];
else
    cpi = [1 3 3];
    snr = [16.5 13.5 13.5];
    elev_ranges{1} = [100 8 6];
    elev_ranges{2} = [100 13 6];
    elev_ranges{3} = [100 22 13 6];
     az_spoil = [1 1 1];
end

for ii = 1:3
    if(~uav )
        [doppler_dwell_time power_dwell_time] = gmtiscript(freqs(ii),snr(ii),cpi(ii),[],[],[],[],[],[],[],[],[],[],[],[], az_spoil(ii),elev_ranges{ii});
    else
        [doppler_dwell_time power_dwell_time] = gmtiscript(freqs(ii),snr(ii),cpi(ii),150,5,[7 100],.5,.2,[],[],[],[],[],[],[],az_spoil(ii),elev_ranges{ii});
    end
    ppt(count) = sum(power_dwell_time);
    ddt(count) = sum(doppler_dwell_time);
    label_text{ii} = sprintf('%3.1f GHz, %3.1f SNR, %d CPI, %d X AZ SPOIL \n',freqs(ii)/1E9,snr(ii),cpi(ii),az_spoil(ii));
    count = count + 1;


end
count = 1;
freqs = [5.5 10 17] * 1E9;
cpi = [1 1 1];
snr = [16.5 16.5 16.5 ];
az_spoil = [1 1 1];
for ii = 1:3
    if(~uav )
        [doppler_dwell_time power_dwell_time] = gmtiscript(freqs(ii),snr(ii),cpi(ii),[],[],[],[],[],[],[],[],[],[],[],[], az_spoil(ii),elev_ranges{ii});
    else
        [doppler_dwell_time power_dwell_time] = gmtiscript(freqs(ii),snr(ii),cpi(ii),150,5,[7 100],.5,.2);
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

legend(label_text{1},label_text{2} ,label_text{3} ,'STAR LEVEL 2 BOUND','POWER BOUND','DOPPLER BOUND','location','north')




hold on
plot(freqs/1E9,att,'g','linewidth',2)
grid on
xlabel('FREQUENCY GHz')
ylabel('TOTAL DWELL TIME')
if(~uav )
    label = sprintf('C:/figures/time_spoiled_elev_beam_isr');
else
    label = sprintf('C:/figures/time_spoiled_beam_uav');
end

print( gcf, '-djpeg', label)
