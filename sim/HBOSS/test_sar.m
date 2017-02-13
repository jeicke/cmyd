%% plot power
% gmtiscript(operating_frequency,SNR,ncpi,Ptransmit,...
%                                                                   radar_height,minmaxrange,width,height,NF,losses,...
%                                                                   T,rcs,azimuth_search,spoil_parameter,...
%                                                                   min_detectable_velocity)
% r3ddw4rf88



width = [2.6 7.3 .5 .6];
height = [.4 .75 .2 .6];
power = [1610 300 165 60];
altitude = [20 20 5 5 ];
range = [25 300;25 300;6 100;6 100];

platform_velocity = [180 200 90 300];
for ii  = 1:length(width)
    
    if(ii<3)
        
        f = [5.5 8 10 14 17] * 1E9;
        cpi = [1 1 1 1 1];
        snr = [25 25 25 25 25];
    else
        f = [5.5 8 10 14 17] * 1E9;
        cpi = [1 1 1 1 1];
        snr = [25 25 25 25 25];
    end
    
    for jj = 1:length(f)
      
        [sar_dwell_time power_dwell_time] = sarscript(f(jj),snr(jj),cpi(jj),power(ii),altitude(ii),range(ii,:),width(ii),height(ii),[],[],[],[],[],[],[],[],[],[],false,platform_velocity(ii) );
        ppt(ii,jj) = sum(power_dwell_time);
        sst(ii,jj) = sum(sar_dwell_time);
        
    end
end
%%
plot(f/1E9,ppt(1,:),'r-x',f/1E9,ppt(2,:),'b-o',f/1E9,ppt(3,:),'g-d',f/1E9,ppt(4,:),'c-p','linewidth',2)
xlabel('FREQUENCY (GHz)')
ylabel('TOTAL TIME FOR SNR')
legend('ISR','JSTARS','UAV','FIGHTER','location','northeast')
grid on
label = sprintf('C:/figures/power_summary_all');
print( gcf, '-djpeg', label )
figure
plot(f/1E9,(sst(1,:)),'r-x',f/1E9,(sst(2,:)),'b-o',f/1E9,(sst(3,:)),'g-d',f/1E9,(sst(4,:)),'c-p','linewidth',2)
xlabel('FREQUENCY (GHz)')
ylabel('TIME FOR 1 M OF SAR RESOLUTION (dBS)')
legend('ISR','JSTARS','UAV','FIGHTER','location','northwest')
grid on
label = sprintf('C:/figures/doppler_summary_all');
print( gcf, '-djpeg', label )