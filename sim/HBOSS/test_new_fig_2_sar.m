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
beams = [32 85 6 7];
range = [25 300;25 300;7 100;7 100];
platform_velocity = [180 200 90 300];
for ii  = 1:4
    
    if(ii<3)
        
        f = [5.5 10 17] * 1E9;
        cpi = [1 1 1];
        snr = [25 25 25];
    else
        f = [5.5 10 17] * 1E9;
        cpi = [1 1 1];
        snr = [25 25 25];
    end
    switch ii
        case 1
            elev_ranges{1} = [ 300.0000   30.1650   25];
            elev_ranges{2} = [300.0000   47.4227   26.1997   25];
            elev_ranges{3} = [300.0000   74.0697   43.2877   30.6287   25];
        case 2
            
            elev_ranges{1} = [ 300.0000   30.1650   25];
            elev_ranges{2} = [300.0000   47.4227   26.1997   25];
            elev_ranges{3} = [300.0000   74.0697   43.2877   30.6287   25];
        case 3
            
            elev_ranges{1} = [100.0000    7.8675    7.0 ];
            elev_ranges{2} = [100.0000   12.6525    7.0];
            elev_ranges{3} = [100.0000   20.4346   11.4874    7.9934    7.0];
        case 4
             elev_ranges{1} = [100.0000    7.8675    7.0 ];
            elev_ranges{2} = [100.0000   12.6525    7.0];
            elev_ranges{3} = [100.0000   20.4346   11.4874    7.9934    7.0];
        otherwise
            disp('Unknown method.')
    end
    
    
    for jj = 1:length(f)
      
        [doppler_dwell_time power_dwell_time az_beams] = sarscript(f(jj),snr(jj),cpi(jj),power(ii),...
                altitude(ii),range(ii,:),width(ii),height(ii),[],[],[],[],[],[],[],[],[],[],[],platform_velocity(ii));
        ppt_1(ii,jj) = sum(power_dwell_time);
        ddt_1(ii,jj) = sum(doppler_dwell_time);
        
        ratio = ceil(sum(doppler_dwell_time)/sum(power_dwell_time));
        ratio2 = ceil(doppler_dwell_time./power_dwell_time);
        if(ratio>beams(ii))
            ratio = beams(ii);
        end
        [doppler_dwell_time power_dwell_time] = sarscript(f(jj),snr(jj),cpi(jj),power(ii),altitude(ii),...
                                                           range(ii,:),width(ii),height(ii),[],[],[],[],[],[],[],[],[],ratio,[],platform_velocity(ii));
        ppt_2(ii,jj) = sum(power_dwell_time);
        ddt_2(ii,jj) = sum(doppler_dwell_time);
        
        ratio = ratio2;
        ratio(find(ratio >beams(ii))) = beams(ii);
        spoil = ratio(floor(length(az_beams)/2):length(az_beams):end);
        
        [doppler_dwell_time power_dwell_time az_beams] = sarscript(f(jj),snr(jj),cpi(jj),power(ii),altitude(ii),...
                                                           range(ii,:),width(ii),height(ii),[],[],[],[],[],[],[],[],[],spoil,[],platform_velocity(ii));
        ppt_3(ii,jj) = sum(power_dwell_time);
        ddt_3(ii,jj) = sum(doppler_dwell_time);
        
        [doppler_dwell_time power_dwell_time az_beams] = sarscript(f(jj),snr(jj),cpi(jj),power(ii),altitude(ii),...
                                                           range(ii,:),width(ii),height(ii),[],[],[],[],[],[],[],[],elev_ranges{jj},[],[],platform_velocity(ii));
        ratio = ceil(doppler_dwell_time./power_dwell_time);
        ratio(find(ratio >beams(ii))) = beams(ii);
        spoil = ratio(floor(length(az_beams)/2):length(az_beams):end);
        [doppler_dwell_time power_dwell_time az_beams] = sarscript(f(jj),snr(jj),cpi(jj),power(ii),altitude(ii),...
                                                           range(ii,:),width(ii),height(ii),[],[],[],[],[],[],[],[],elev_ranges{jj},spoil,[],platform_velocity(ii));
        ppt_4(ii,jj) = sum(power_dwell_time);
        ddt_4(ii,jj) = sum(doppler_dwell_time);
                                                       
        [doppler_dwell_time power_dwell_time az_beams] = sarscript(f(jj),16.5,1,power(ii),altitude(ii),...
                                                           range(ii,:),width(ii),height(ii),[],[],[],[],[],[],[],[],elev_ranges{jj},[],[],platform_velocity(ii));
        ratio = ceil(doppler_dwell_time./power_dwell_time);
        ratio(find(ratio >beams(ii))) = beams(ii);
        spoil = ratio(floor(length(az_beams)/2):length(az_beams):end);
        [doppler_dwell_time power_dwell_time az_beams] = sarscript(f(jj),16.5,1,power(ii),altitude(ii),...
                                                           range(ii,:),width(ii),height(ii),[],[],[],[],[],[],[],[],elev_ranges{jj},spoil,[],platform_velocity(ii));
        ppt_5(ii,jj) = sum(power_dwell_time);
        ddt_5(ii,jj) = sum(doppler_dwell_time);                                           
    end
end
%%
x = [1:4];
for index = 1:4;
single_beam = max([ddt_1(index,:); ppt_1(index,:)]);
multi_beam_az_avg = max([ddt_2(index,:); ppt_2(index,:)]);
multi_beam_az = max([ddt_3(index,:); ppt_3(index,:)]);
multi_beam_az_el = max([ddt_4(index,:); ppt_4(index,:)]);
multi_beam_az_el_star = max([ddt_5(index,:); ppt_5(index,:)]);
figure;
p1 = [single_beam(1) multi_beam_az(1) multi_beam_az_el(1) multi_beam_az_el_star(1)];
p2 = [single_beam(2) multi_beam_az(2) multi_beam_az_el(2) multi_beam_az_el_star(2)];
p3 = [single_beam(3) multi_beam_az(3) multi_beam_az_el(3) multi_beam_az_el_star(3)];
plot(x,p1,'gx-',x,p2,'ro-',x,p3,'cd-','linewidth',2);
switch index
    case 1;
        ylabel('ISR TIME (SEC)')
    case 2;
        ylabel('JSTARS TIME (SEC)')
    case 3;
        ylabel('UAV TIME (SEC)')
    case 4;
        ylabel('FIGHTER TIME (SEC)')
end
        
   
legend('C','X','KU');

set(gca,'xtick',[1 2 3 4 ]);
set(gca,'XTickLabel',['SINGLE BEAM';
                      'MULT AZ    ';
                      'MULT AZ& EL';
                      'LVL 2 STAR '])
grid;             
end