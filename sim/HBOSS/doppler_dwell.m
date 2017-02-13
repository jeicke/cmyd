clear all;
width = 2.4/2; %in m
height = .4/2; % in m
res = 1;
min_doppler_resolution = 1; % in m/s radial velocity
%% radar equation parameters

NF=3;
losses=5;
Ptav=1100;
lamd=0.03;
sigma=0;
SNR=25;
k  = 1.3806503*10^-23;
T = 290;
%% operating parameters
freqs = [5 10 19] * 1E9;
index2 = 1;
for operating_frequency  = freqs;

c = 299792458; % speed of light in m/s
index = 1;
rheight = 5:.1:15;
for radar_height = rheight; % in KM
min_range = 25; % in km
max_range = 300; %in km
azimuth_search = [-60 60];
% calculate beamsize in degrees
wavelength = c/operating_frequency;
beamfactor=17/20; %factor to convert to 3dB beamwidth

azbw = asind(2 * wavelength/width) * beamfactor*.5;
elbw = asind(2 * wavelength/height) * beamfactor*.5;
% create example array 1/2 wavelength spaced
count = 1;
for w = [-width/2:wavelength/2:width/2]
    for h = [-height/2:wavelength/2:height/2]
        geometry(count,1) =w;
        geometry(count,2) = 0;
        geometry(count,3)= h;
        count = count +1;
    end
end
% calculate plane to ground angles and ground to plane angels

count = 1;
clear ele;
for range = min_range:max_range
    [elevation,grazing]=sargmtiangles(radar_height,range);
    graze(count) = grazing;
    elev(count) = elevation;
    count = count + 1;
end
range = min_range:max_range;
% create beam pointing centers
count = 1;

beams =(min(elev)):elbw/2:(max(elev));

for beamcenter_el = beams 

  for beamcenter_az = [azimuth_search(1):azbw/2:azimuth_search(2)]
      beam_az(count) = beamcenter_az;
      beam_el(count) = beamcenter_el;

      beam_ranges(count)=sargmtiangles2(radar_height,beam_el(count));

      beam_x(count) = beam_ranges(count) * sind(beam_az(count));
      beam_y(count) =  beam_ranges(count) * cosd(beam_az(count));
      count = count + 1;
  end
end
% compute dwell time for each beam based on Doppler Resolution

doppler_shift(1) =operating_frequency- operating_frequency * (1-min_doppler_resolution/c);
% using rayleigh frequency resolution
doppler_dwell_time = 1/doppler_shift(1);
doppler_dwell_time = wavelength/min_doppler_resolution * 1./sind(90-max(abs(elev)));
doppler_dwell_time_beam = wavelength/min_doppler_resolution * 1./sind(90-(abs(beam_el)));
total_doppler_dwell_time(index,index2) = sum(doppler_dwell_time_beam);
mean_doppler_dwell_time(index,index2) = mean(doppler_dwell_time_beam);
index = index + 1;
end
index2 = index2+1
end
%%
plot(rheight,total_doppler_dwell_time(:,1),'r',rheight,total_doppler_dwell_time(:,2),'b',rheight,total_doppler_dwell_time(:,3),'g','linewidth',2)

legend('5 GHz', '10 GHz','19 GHz','Location','Best');
xlabel('RADAR ELEVATION (KM)')
ylabel('TOTAL DWELL FOR 1 M/S DOPPLER RESOLUTION (SEC)')
grid
title('TIME NEEDED FOR DOPPLER RESOLUTION SPOILED BEAMS 100-300 KM')
