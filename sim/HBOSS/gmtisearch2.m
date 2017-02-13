%% array size parameters
clear all;
indexs = 1;
tf =  [ 5E9  8E9 10E9 12E9 15E9 17E9];
for operating_frequency = tf
fprintf(1,'OPERATING FREQUENCY %d\n',operating_frequency/1E9);
% width = .5; %in m
% height = .2; % in m
%% physical parameters
c = 299792458; % speed of light in m/s
k  = 1.3806503*10^-23; % Boltzmann's constant
%% radar equation parameters


width = 2.4; %in m
height = .4; % in m
SNR=10.5;
res = 1;
NF=3;
losses=5;
Ptransmit=1500;
T = 290;
rcs = 1;

radar_height = 20; % in KM
min_range = 25; % in km
max_range = 300; %in km
azimuth_search = [-60 60];
spoil_parameter =1;
%% calculate beamsize in degrees receive
wavelength = c/operating_frequency;
beamfactor=17/20; %factor to convert to 3dB beamwidth

azbw = asind(2 * wavelength/width) * beamfactor;
elbw = asind(2 * wavelength/height) * beamfactor;
%% calculate beamsize in degrees transmit
wavelength = c/operating_frequency;
beamfactor=17/20; %factor to convert to 3dB beamwidth

azbw_transmit = asind(2 * wavelength/(width * spoil_parameter)) * .5 * beamfactor ;
elbw_transmit = asind(2 * wavelength/(height * spoil_parameter)) * .5 * beamfactor ;
%% create example array 1/2 wavelength spaced receive
count = 1;
wl = c/operating_frequency;
clear geometry;
for w = [-width/2:wl/2:width/2]
    for h = [-height/2:wl/2:height/2]
        geometry(count,1) =w;
        geometry(count,2) = 0;
        geometry(count,3)= h;
        count = count +1;
    end
end
%% create example array 1/2 wavelength spaced transmit
count = 1;
wl =  c/operating_frequency;
clear geometry_transmit;
for w = [-(width * spoil_parameter)/2:wl/2:(width * spoil_parameter)/2]
    for h = [-(height * spoil_parameter)/2:wl/2:(height * spoil_parameter)/2]
        geometry_transmit(count,1) =w;
        geometry_transmit(count,2) = 0;
        geometry_transmit(count,3)= h;
        count = count +1;
    end
end
%% calculate gain of array
Gr=41.5;%+20*log10(3);
Gt=41.5;
clear HGA;
HGA.theta_source =0;
HGA.phi_source = 0;
HGA.SLx = geometry;
resolution = 1;
%[theory_matrix az de HGA] = calcBeamPattern(operating_frequency,[0:resolution :360],[-90:resolution:89],[],HGA);
%[DI PSG] = calculatePSG(resolution,theory_matrix,az,de,ones(size(de)))

PSG = 10*log10(4 * pi * width*height/wavelength^2);
PSG_transmit = 10*log10(4 * pi * width*height *  spoil_parameter^2/wavelength^2);
range = min_range:max_range;
% %% calculate plane to ground angles and ground to plane angels
% 
% count = 1;
% clear ele;
% %go max to min
% [elev,graze]=groundsurv(radar_height,min_range,max_range, sind(elbw));
% % for range = min_range:max_range
% %     [elevation,grazing]=sargmtiangles(radar_height,range);
% %     graze(count) = grazing;
% %     elev(count) = elevation;
% %     count = count + 1;
% % end

%% create transmit beam pointing centers
count = 1;
[elev,graze]=groundsurv(radar_height,min_range,max_range, sind(elbw_transmit));
beams =elev;%(min(elev)):elbw_transmit/2:(max(elev));

for beamcenter_el = beams 

  for beamcenter_az = [azimuth_search(1):azbw_transmit/2:azimuth_search(2)]
      beam_az_transmit(count) = beamcenter_az;
      beam_el_transmit(count) = beamcenter_el;

      beam_ranges_transmit(count)=sargmtiangles2(radar_height,beam_el_transmit(count));
      beam_max_range_transmit(count) = sargmtiangles2(radar_height,beam_el_transmit(count) + elbw_transmit/2);
      if( abs( beam_max_range_transmit(count))~= beam_max_range_transmit(count))
           beam_max_range_transmit(count) = max_range;
      end
      beam_x_transmit(count) = beam_ranges_transmit(count) * sind(beam_az_transmit(count));
      beam_y_transmit(count) =  beam_ranges_transmit(count) * cosd(beam_az_transmit(count));
      count = count + 1;
  end
end
%% create receive beam pointing centers
count = 1;
[elev,graze]=groundsurv(radar_height,min_range,max_range, sind(elbw));
beams =elev;

for beamcenter_el = beams 

  for beamcenter_az = [azimuth_search(1):azbw/2:azimuth_search(2)]
      beam_az(count) = beamcenter_az;
      beam_el(count) = beamcenter_el;

      beam_ranges(count)=sargmtiangles2(radar_height,beam_el(count));
      beam_max_range(count) = sargmtiangles2(radar_height,beam_el(count) + elbw/2);
      if( abs(beam_max_range(count))~=beam_max_range(count))
          beam_max_range(count) = max_range;
      end
      beam_x(count) = beam_ranges(count) * sind(beam_az(count));
      beam_y(count) =  beam_ranges(count) * cosd(beam_az(count));
      count = count + 1;
  end
end
%% compute dwell time for each beam based on Doppler Resolution
min_doppler_resolution = 1; % in m/s radial velocity
doppler_shift(1) =operating_frequency- operating_frequency * (1-min_doppler_resolution/c);
% using rayleigh frequency resolution
doppler_dwell_time = 1/doppler_shift(1);
doppler_dwell_time = wavelength/(2*min_doppler_resolution) * 1./sind(90-max(abs(elev)));
doppler_dwell_time_beam = wavelength/(2*min_doppler_resolution) * 1./sind(90-(abs(beam_el_transmit)));
total_doppler_dwell_time = sum(doppler_dwell_time_beam);
fprintf(1,'TOTAL DWELL TIME NEEDED FOR DOPPLER %3.1f\n',sum(total_doppler_dwell_time));
ddt(indexs ) = sum(total_doppler_dwell_time);

max_beams = 10/doppler_dwell_time;

% %% beam power
% for ii = 1:length(beam_max_range)
%     Ptav(ii) = (4 * pi)^3 * (beam_max_range(ii) * 1000).^4 * 10^(NF/10) * k * T * 10^(losses/10) * 10^(SNR/10)./...
%                     (10.^(PSG_transmit/10)  * doppler_dwell_time_beam(ii) *10^(PSG/10).* cosd(beam_az(ii)).^2 .* cosd(beam_el(ii)).^2);
% end
%% compute dwell time based on power


power_dwell_time = (4 * pi)^2 * (beam_max_range_transmit * 1000).^4 * 10^(NF/10) * k * T * 10^(losses/10) * 10^(SNR/10)./...
                    (10.^(PSG/10) * Ptransmit *   rcs * width * height .* cosd(beam_az_transmit).^2 .* cosd(beam_el_transmit).^2);
                
fprintf(1,'TOTAL DWELL TIME NEEDED FOR POWER %3.1f\n',sum(power_dwell_time));

ppt(indexs ) = sum(power_dwell_time);
max_beams = 10/sum(power_dwell_time);
clear power_dwell_time;
indexs  = indexs +1;
fprintf(1,'\n');


%%
% clear X;
% clear Y;
% clear Ptav;
% clear scallop_transmit;
% clear scallop_receive;
% clear doppler_dwell_time;
% clear power_dwell_time;
% count1 = 1;
% count2 =  1;
% res = 1;
% trans_hga.SLx = geometry_transmit;
% trans_hga.usestored = false;
% rec_hga.SLx = geometry;
% rec_hga.usestored = false;
% use_scallop_loss = false;
% clear power_dwell_time;
% for range = min_range:res:max_range
%     for az = [azimuth_search(1):res:azimuth_search(2)]
%         range_all(count2,count1) = range;
%         X(count2,count1) = range * sind(az);
%         Y(count2,count1) = range * cosd(az);
%         [el,grazing]=sargmtiangles(radar_height,range);
% 
%         % compute scalloping loss transmit and receive
% 
%         %find closes transmit beam
%         range_from_beam = sqrt( (X(count2,count1) - beam_x_transmit).^2 +  (Y(count2,count1) - beam_y_transmit).^2 );
% 
%         [j closest_beam]  = min(range_from_beam);
%         closest_beam = closest_beam(1);
%         trans_hga.theta_source =beam_az_transmit(closest_beam);
%         trans_hga.phi_source = beam_el_transmit(closest_beam);
% 
%         [scallop_transmit(count2,count1)] = double(calcBeamPattern(operating_frequency,az,el,[],trans_hga));
%         range_from_beam = sqrt( (X(count2,count1) - beam_x).^2 +  (Y(count2,count1) - beam_y).^2 );
% 
%         [j closest_beam]  = min(range_from_beam);
%         closest_beam = closest_beam(1);
%         rec_hga.theta_source =beam_az(closest_beam);
%         rec_hga.phi_source = beam_el(closest_beam);
%         [scallop_receive(count2,count1)] = double(calcBeamPattern(operating_frequency,az,el,[],rec_hga));
% 
% 
%         power_dwell_time(count2,count1) = (4 * pi)^2 * (range  * 1000).^4 * ...
%             10^(NF/10) * k * T  * 10^((losses-2)/10) * 10^(SNR/10)./...
%             (Ptransmit * rcs * 10.^(PSG_transmit/10) *...
%             (scallop_transmit(count2,count1)) * width * height * scallop_receive(count2,count1) ...
%             .* cosd(az).^2 .* cosd(el).^2 );
%         power_dwell_time_noscallop(count2,count1) = (4 * pi)^2 * (range  * 1000).^4 * ...
%             10^(NF/10) * k * T * 10^(losses/10) * 10^(SNR/10)./...
%             (Ptransmit * rcs * width * height *10.^(PSG_transmit/10) ...
%             .* cosd(az).^2 .* cosd(el).^2);
%         
%         
%         doppler_dwell_time(count2,count1)  = wavelength/min_doppler_resolution * 1./sind(90-max(abs(el)));
%         count1 = count1 + 1;
%     end
%     count1 = 1;
%     count2 = count2 + 1;
% end
% %%
% figure
% pcolor(X,Y,10*log10(scallop_transmit));
% shading flat;
% axis equal
% axis tight
% xlabel('Y KM')
% ylabel('X KM')
% caxis([-5 0]);
% colorbar('location','southoutside')
% hold on
% plot(beam_x_transmit,beam_y_transmit,'xw')
% title('TRANSMIT SCALLOP')
% figure
% pcolor(X,Y,10*log10(scallop_receive));
% shading flat;
% axis equal
% axis tight
% xlabel('Y KM')
% ylabel('X KM')
% caxis([-5 0]);
% colorbar('location','southoutside')
% hold on;
% plot(beam_x,beam_y,'dw')
% title('RECEIVE SCALLOP')
% 
% 
% figure
% pcolor(X,Y,10*log10(power_dwell_time));
% shading flat;
% axis equal
% axis tight
% xlabel('Y KM')
% ylabel('X KM')
% colorbar('location','southoutside')
% title('POWER DWELL TIME dBs')
% figure
% pcolor(X,Y,10*log10(power_dwell_time_noscallop));
% shading flat;
% axis equal
% axis tight
% xlabel('Y KM')
% ylabel('X KM')
% colorbar('location','southoutside')
% title('POWER DWELL TIME dBs NO SCALLOP')
% figure
% pcolor(X,Y,10*log10(doppler_dwell_time));
% shading flat;
% axis equal
% axis tight
% xlabel('Y KM')
% ylabel('X KM')
% colorbar('location','southoutside')
% title('DOPPLER DWELL TIME dBs')
% figure
% 
% 
% pcolor(X,Y,10*log10(power_dwell_time./doppler_dwell_time));
% shading flat;
% axis equal
% axis tight
% xlabel('Y KM')
% ylabel('X KM')
% h = colorbar('location','southoutside');
% hold on;
% %plot(beam_x_transmit,beam_y_transmit,'xw','markersize',.5)
% %plot(beam_x,beam_y,'ow','markersize',1)
% title('RATIO OF POWER AND DOPPLER DWELL TIMES')
% xlabel(h,'(TIME NEEDED FOR POWER) / (TIME NEEDED FOR DOPPLER) dB');
% label = sprintf('C:/figures/ratio_ns_%dghz_uav',operating_frequency/1e9);
% caxis([-40 15])
% print( gcf, '-djpeg', label )
% 
% figure
% pcolor(X,Y,10*log10(power_dwell_time_noscallop./doppler_dwell_time));
% shading flat;
% axis equal
% axis tight
% xlabel('Y KM')
% ylabel('X KM')
% h = colorbar('location','southoutside');
% title('RATIO OF POWER AND DOPPLER DWELL TIMES W/O SCALLOP')
% xlabel(h,'(TIME NEEDED FOR POWER) / (TIME NEEDED FOR DOPPLER) dB');
% caxis([-40 15])
% label = sprintf('C:/figures/ratio_sc_%dghz_uav',operating_frequency/1e9);
% print( gcf, '-djpeg', label )
end
%%
plot(tf/1E9,ppt,'r-x','linewidth',2)
xlabel('FREQUENCY GHz')
ylabel('TOTAL DWELL TIME FOR 10.5dB of SNR')
grid on
label = sprintf('C:/figures/power_summary_uav');
print( gcf, '-djpeg', label )
figure
plot(tf/1E9,ddt,'r-x','linewidth',2)
xlabel('FREQUENCY GHz')
ylabel('TOTAL DWELL TIME FOR 1 M/S OF DOPPLER RESOLUTION')
grid on
label = sprintf('C:/figures/doppler_summary_uav');
print( gcf, '-djpeg', label )
%% create beam in grid
% clear HGA
% for ii = 1:length(beam_ranges)
%     range = beam_ranges(ii);
% 
%     [elevation,grazing]=SARGMTIangles(radar_height,range );
%     HGA.theta_source =beam_az(ii);
%     HGA.phi_source = elevation;
%     HGA.SLx = geometry;
%     HGA.usestored = true;
%     [theory_matrix az de HGA] = calcBeamPattern(operating_frequency,[azimuth_search(1):res:azimuth_search(2)],elev,[],HGA);
% 
%     count2 = 1;
%     count1 = 1;
%     clear X;
%     clear Y;
%     for range = min_range:max_range
%         for az = [azimuth_search(1):res:azimuth_search(2)]
%             range_all(count2,count1) = range;
%             X(count2,count1) = range * sind(az);
%             Y(count2,count1) = range * cosd(az);
%             count1 = count1 + 1;
%         end
%         count1 = 1;
%         count2 = count2 + 1;
%     end
%     hold off
% 
%     Ptav = (4 * pi)^3 * (range_all * 1000).^4 * 10^(NF/10) * k * T * 10^(losses/10) * 10^(SNR/10)./...
%         (10.^(PSG/10)  * power_dwell_time(1) *10^(PSG/10))
% 
% 
%     Ptav(:,:,ii) = Ptav;
%     pcolor(X,Y,Ptav);
%     shading flat;
%     axis equal
%     axis tight
%     xlabel('Y KM')
%     ylabel('X KM')
%     caxis([max(max(Ptav))-3 max(max(Ptav)) ]);
%     hold on
%     plot(beam_x,beam_y,'xw')
%     colorbar;
%     drawnow
% end
% %%
% snr_mean = 10*log10(sum(10.^(snr_all/10),3));
% pcolor(X,Y,snr_mean);
% shading flat;
% axis equal
% axis tight
% xlabel('Y KM')
% ylabel('X KM')
% caxis([max(max(snr_mean))-100 max(max(snr_mean)) ]);
% hold on
% plot(beam_x,beam_y,'xw')
% colorbar;


 