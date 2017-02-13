function [beam_transmit beam_receive doppler_dwell_time_beam ] = gmtigrid(operating_frequency,SNR,Ptransmit,...
                                                                          radar_height,minmaxrange,width,height,...
                                                                          NF,losses,T,rcs,azimuth_search,spoil_parameter,...
                                                                          min_detectable_velocity,sim_recv_beams,...
                                                                          elev_beam_ranges,spoil)

if(~exist('operating_frequency','var')||isempty(operating_frequency))
    operating_frequency = 5E9; %center frequency in Hz
end
if(~exist('SNR','var')||isempty(SNR))
    SNR = 10.5; %SNR in dB
end

if(~exist('Ptransmit','var')||isempty(Ptransmit))
    Ptransmit = 1500; % transmitted power in Watts
end
if(~exist('radar_height','var')||isempty(radar_height))
    radar_height = 20; % radar height in km
end
if(~exist('minmaxrange','var')||isempty(minmaxrange))
    min_range = 25; % radar height in km
    max_range = 300;
else
    min_range = minmaxrange(1); % radar height in km
    max_range = minmaxrange(2);
end
if(~exist('width','var')||isempty(width))
        width = 2.4; % width of array in meters
end
if(~exist('height','var')||isempty(height))
        height = .4; % height of array in meters
end
if(~exist('NF','var')||isempty(NF))
        NF = 3; %receiver noise figure in dB
end
if(~exist('losses','var')||isempty(losses))
    losses = 5; % losses in dB
end
if(~exist('rcs','var')||isempty(rcs))
    rcs = 1; % radar cross section in m^2
end
if(~exist('T','var')||isempty(T))
    T = 290; % receiver temperature in Kelvin
end
if(~exist('azimuth_search','var')||isempty(azimuth_search))
    azimuth_search = [-60 60]; % azimuthal scan window in start and end degrees
end
if(~exist('spoil_parameter','var')||isempty(spoil_parameter))
    spoil_parameter = 1; % transmit spoil in fractional amount
end
if(~exist('min_detectable_velocity','var')||isempty(min_detectable_velocity))
    min_detectable_velocity = 1; % min detectable velocity in m/s
end
if(~exist('sim_recv_beams','var')||isempty(sim_recv_beams))
    sim_recv_beams = 1;
end

if(~exist('elev_beam_ranges','var')||isempty(elev_beam_ranges))
    elev_beam_ranges = [];
end
if(~exist('spoil','var')||isempty(spoil))
    spoil = [];
end
% physical parameters
c = 299792458; % speed of light in m/s
k  = 1.3806503*10^-23; % Boltzmann's constant

% calculate beamsize in degrees
wavelength = c/operating_frequency;
beamfactor=17/20; %factor to convert to 3dB beamwidth

azbw_transmit = asind(wavelength/width) * beamfactor ;
elbw_transmit = asind(wavelength/height) * beamfactor;

azbw = asind(wavelength/width) * beamfactor;
elbw = asind(wavelength/height) * beamfactor *.5;


% create transmit beam pointing centers

if(isempty(elev_beam_ranges))
    [elev,graze,sranges]=groundsurv(radar_height,min_range,max_range, sind(elbw_transmit));
else
    [elev,graze]=SARGMTIangles(radar_height,elev_beam_ranges);
end

beams =elev;%(min(elev)):elbw_transmit/2:(max(elev));
%
if(isempty(spoil))
    % calculate beamsize in degrees transmit
    if(sim_recv_beams>1)
        width2 = wavelength/sin( asin(wavelength/width)*beamfactor*sim_recv_beams);
    else
        width2 = width;
    end

    azbw_transmit = asind(wavelength/(width * spoil_parameter)) * beamfactor  * sim_recv_beams;
    elbw_transmit = asind(wavelength/(height * spoil_parameter))  * beamfactor * .5;


    % calculate gain of array
    PSG = 10*log10(4 * pi * width*height/wavelength^2);
    PSG_transmit = 10*log10(4 * pi * width2*height/wavelength^2);
    range = min_range:max_range;
    count = 1;
    az_beams = [azimuth_search(1):azbw_transmit:(azimuth_search(2)+azbw_transmit/2)];
    count2 = 1;
    for beamcenter_el = beams

        for beamcenter_az = az_beams
            beam_az_transmit(count) = beamcenter_az;
            beam_el_transmit(count) = beamcenter_el- elbw_transmit/2;

            beam_ranges_transmit(count)=sargmtiangles2(radar_height,beam_el_transmit(count));
            beam_max_range_transmit(count) = sargmtiangles2(radar_height,beam_el_transmit(count) + elbw_transmit/2);
            beam_min_range_transmit(count) = sargmtiangles2(radar_height,beam_el_transmit(count) - elbw_transmit/2);
            if( abs( beam_max_range_transmit(count))~= beam_max_range_transmit(count))
                beam_max_range_transmit(count) = max_range;
            end
            if(beam_max_range_transmit(count)>max_range)
                beam_max_range_transmit(count) = max_range;
            end
            if(beam_max_range_transmit(count)<min_range)
                beam_min_range_transmit(count) = min_range;
            end
            PSG_transmit(count) = 10*log10(4 * pi * width2*height/wavelength^2);
            beam_x_transmit(count) = beam_ranges_transmit(count) * sind(beam_az_transmit(count));
            beam_y_transmit(count) =  beam_ranges_transmit(count) * cosd(beam_az_transmit(count));
           
            count = count + 1;
        end
        count2 = count2  + 1;
    end
else
    width2 = wavelength./sin( asin(wavelength/width)*beamfactor.*spoil);
    azbw_transmit = asind(wavelength/(width )) * beamfactor  * spoil;
    count = 1;
    count2 = 1;
    
    for beamcenter_el = beams
        if(length(azbw_transmit)<2)
            az_beams = [azimuth_search(1):azbw_transmit:(azimuth_search(2)+azbw_transmit/2)];
        else
            az_beams = [azimuth_search(1):azbw_transmit(count2):(azimuth_search(2)+azbw_transmit(count2)/2)];
        end
        for beamcenter_az = az_beams
            beam_az_transmit(count) = beamcenter_az;
            beam_el_transmit(count) = beamcenter_el;

            beam_ranges_transmit(count)=sargmtiangles2(radar_height,beam_el_transmit(count));
            beam_max_range_transmit(count) = sargmtiangles2(radar_height,beam_el_transmit(count) + elbw_transmit/2);
            beam_min_range_transmit(count) = sargmtiangles2(radar_height,beam_el_transmit(count) - elbw_transmit/2);
            if(abs(beam_max_range_transmit(count))~= beam_max_range_transmit(count))
                beam_max_range_transmit(count) = max_range;
            end
            if(beam_max_range_transmit(count)>max_range)
                beam_max_range_transmit(count) = max_range;
            end
            if(~isreal(beam_max_range_transmit(count)))
                beam_max_range_transmit(count) = max_range;
            end
            if(beam_max_range_transmit(count)>min_range)
                beam_min_range_transmit(count) = min_range;
            end
            beam_x_transmit(count) = beam_ranges_transmit(count) * sind(beam_az_transmit(count));
            beam_y_transmit(count) =  beam_ranges_transmit(count) * cosd(beam_az_transmit(count));
            if(length(azbw_transmit)<2)
                PSG_transmit(count) = 10*log10(4 * pi * width2*height/wavelength^2);
            else
                PSG_transmit(count) = 10*log10(4 * pi * width2(count2)*height/wavelength^2);
            end
            if(length(spoil)>1)
             recv_beams(count) = spoil(count2);
            else
                recv_beams(count) =  1;
            end
            count = count + 1;
        end
        count2  = count2 + 1;
    end
    
end

% compute dwell time for each beam based on Doppler Resolution
min_doppler_resolution = min_detectable_velocity; % in m/s radial velocity

% using rayleigh frequency resolution
doppler_dwell_time = wavelength/(2*min_doppler_resolution) * 1./sind(90-max(abs(elev)));
doppler_dwell_time_beam = wavelength/(2*min_doppler_resolution) * 1./sind(90-(abs(beam_el_transmit)));

% create receive beams and associate with closest transmit beams
[elev,graze,sranges]=groundsurv(radar_height,min_range,max_range, sind(elbw));
az_beams = [azimuth_search(1):azbw:(azimuth_search(2)+azbw/2)];
count = 1;
for beamcenter_el = beams

    for beamcenter_az = az_beams
        beam_az(count) = beamcenter_az;
        beam_el(count) = beamcenter_el- elbw/2;

        beam_ranges(count)=sargmtiangles2(radar_height,beam_el(count));
        beam_min_range(count) = sargmtiangles2(radar_height,beam_el(count) - elbw/2);
        beam_max_range(count) = sargmtiangles2(radar_height,beam_el(count) + elbw/2);
        if( abs( beam_max_range(count))~= beam_max_range(count))
            beam_max_range_transmit(count) = max_range;
        end
        if(~isreal(beam_max_range(count)))
                beam_max_range(count) = max_range;
            end
        beam_x(count) = beam_ranges(count) * sind(beam_az(count));
        beam_y(count) =  beam_ranges(count) * cosd(beam_az(count));
        range_from_beam = sqrt( (beam_x(count) - beam_x_transmit).^2 +  (beam_y(count) - beam_y_transmit).^2 );
        [j closest_beam]  = min(range_from_beam);
        closest_beam = closest_beam(1);
        beam_closest(count) = closest_beam;
        PSG(count) = 10*log10(4 * pi * width*height/wavelength^2);
        count = count + 1;

    end
end
beam_transmit = [beam_az_transmit; beam_el_transmit;beam_ranges_transmit;beam_max_range_transmit;beam_min_range_transmit;PSG_transmit;doppler_dwell_time_beam ];

beam_receive = [beam_az; beam_el;beam_ranges;beam_max_range;beam_min_range;PSG];