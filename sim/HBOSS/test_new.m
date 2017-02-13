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


for ii  = 1:length(width)
    
    if(ii<3)
        
        f = [5.5 8 10 14 17] * 1E9;
        cpi = [3 4 5 6 7];
        snr = [16.5 16.5 16.5 16.5 16.5];
    else
        f = [5.5 8 10 14 17] * 1E9;
        cpi = [1 2 2 3 3];
        snr = [16.5 16.5 16.5 16.5 16.5];
    end
    
    for jj = 1:length(f)
      
        [doppler_dwell_time power_dwell_time] = gmtiscript(f(jj),snr(jj),cpi(jj),power(ii),altitude(ii),range(ii,:),width(ii),height(ii),[],[],[],[],[],[],[],1);
        ppt(ii,jj) = sum(power_dwell_time);
        ddt(ii,jj) = sum(doppler_dwell_time);
        
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
plot(f/1E9,10*log10(ddt(1,:)),'r-x',f/1E9,10*log10(ddt(2,:)),'b-o',f/1E9,10*log10(ddt(3,:)),'g-d',f/1E9,10*log10(ddt(4,:)),'c-p','linewidth',2)
xlabel('FREQUENCY (GHz)')
ylabel('TIME FOR 1 M/S OF VELOCITY RESOLUTION (dBS)')
legend('ISR','JSTARS','UAV','FIGHTER','location','northwest')
grid on
label = sprintf('C:/figures/doppler_summary_all');
print( gcf, '-djpeg', label )