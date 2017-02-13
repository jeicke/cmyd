%% array size parameters
clear all;
width = 2.4; %in m
height = .4; % in m
res = 1;
%% radar equation parameters

NF=3;
losses=5;
Ptav=100;
lamd=0.03;
sigma=0;
SNR=25;
k  = 1.3806503*10^-23;
T = 290;
%% operating parameters
operating_frequency = 5E9;  % in Hz
c = 299792458; % speed of light in m/s
radar_height = 5; % in KM
min_range = 24; % in km
max_range = 100; %in km
azimuth_search = [-60 60];
%% calculate beamsize in degrees
wavelength = c/operating_frequency;
beamfactor=17/20; %factor to convert to 3dB beamwidth

azbw = asind(2 * wavelength/width) * beamfactor;
elbw = asind(2 * wavelength/height) * beamfactor;
%% set if looking up or down
% if true we need to use grazing angles to steer beams, otherwise 
looking_down = true;


%% create example array 1/2 wavelength spaced
count = 1;
for w = [-width/2:wavelength/2:width/2]
    for h = [-height/2:wavelength/2:height/2]
        geometry(count,1) =w;
        geometry(count,2) = 0;
        geometry(count,3)= h;
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
% better way to calculate gain?
%[theory_matrix az de HGA] = calcBeamPattern(operating_frequency,[0:resolution :360],[-90:resolution:89],[],HGA);
%[DI PSG] = calculatePSG(resolution,theory_matrix,az,de,ones(size(de)))

PSG = 10*log10(4 * pi * width*height/wavelength^2)
%% calculate plane to ground angles and ground to plane angels

count = 1;
clear ele;
for range = min_range:max_range
    [elevation,grazing]=sargmtiangles(radar_height,range);
    graze(count) = grazing;
    elev(count) = elevation;
    count = count + 1;
end
range = min_range:max_range;
%% create beam pointing centers
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
%% compute dwell time for each beam based on Doppler Resolution
min_doppler_resolution = 1; % in m/s radial velocity
doppler_shift(1) =operating_frequency- operating_frequency * (1-min_doppler_resolution/c)
% using rayleigh frequency resolution
doppler_dwell_time = 1/doppler_shift(1);
doppler_dwell_time = wavelength/min_doppler_resolution * 1./sind(90-max(abs(elev)))
doppler_dwell_time_beam = wavelength/min_doppler_resolution * 1./sind(90-(abs(beam_el)))
total_doppler_dwell_time = sum(doppler_dwell_time_beam);
length(beam_ranges)
max_beams = 10/doppler_dwell_time
%% compute dwell time based on power

power_dwell_time = (4 * pi)^3 * (beam_ranges * 1000).^4 * 10^(NF/10) * k * T * 10^(losses/10) * 10^(SNR/10)./...
                    (10.^(PSG/10) * Ptav * 10^(PSG/10));
max_beams = 10/sum(power_dwell_time)
%% create beam in grid
clear HGA
power_dwell_time = 10;
for ii = 1:length(beam_ranges)
    range = beam_ranges(ii);

    [elevation,grazing]=SARGMTIangles(radar_height,range );
    HGA.theta_source =beam_az(ii);
    HGA.phi_source = elevation;
    HGA.SLx = geometry;
    HGA.usestored = true;
    [theory_matrix az de HGA] = calcBeamPattern(operating_frequency,[azimuth_search(1):res:azimuth_search(2)],elev,[],HGA);

    count2 = 1;
    count1 = 1;
    clear X;
    clear Y;
    for range = min_range:max_range
        for az = [azimuth_search(1):res:azimuth_search(2)]
            range_all(count2,count1) = range;
            X(count2,count1) = range * sind(az);
            Y(count2,count1) = range * cosd(az);
            count1 = count1 + 1;
        end
        count1 = 1;
        count2 = count2 + 1;
    end
    hold off
    
    Ptav = (4 * pi)^3 * (range_all * 1000).^4 * 10^(NF/10) * k * T * 10^(losses/10) * 10^(SNR/10)./...
                    (10.^(PSG/10)  * power_dwell_time(1) *10^(PSG/10))
    
    
    Ptav(:,:,ii) = Ptav;
    pcolor(X,Y,Ptav);
    shading flat;
    axis equal
    axis tight
    xlabel('Y KM')
    ylabel('X KM')
    caxis([max(max(Ptav))-3 max(max(Ptav)) ]);
    hold on
    plot(beam_x,beam_y,'xw')
    colorbar;
    drawnow
end
%% 
snr_mean = 10*log10(sum(10.^(snr_all/10),3));
 pcolor(X,Y,snr_mean);
shading flat;
    axis equal
    axis tight
    xlabel('Y KM')
    ylabel('X KM')
    caxis([max(max(snr_mean))-100 max(max(snr_mean)) ]);
    hold on
    plot(beam_x,beam_y,'xw')
    colorbar;
 

 