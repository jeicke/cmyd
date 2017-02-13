% array size parameters
function [doppler_dwell_time_beam power_dwell_time az_beams power_dwell_time_band doppler_dwell_time_band   X  Y ] = ...
                                                                  gmtiscript(operating_frequency,SNR,ncpi,Ptransmit,...
                                                                  radar_height,minmaxrange,width,height,NF,losses,...
                                                                  T,rcs,azimuth_search,spoil_parameter,...
                                                                  min_detectable_velocity,sim_recv_beams,elev_beam_ranges,spoil,...
                                                                  compute_bandgraph)
% set parameters
doppler_dwell_time_beam= [];
power_dwell_time =[];
if(~exist('operating_frequency','var')||isempty(operating_frequency))
    operating_frequency = 5E9; %center frequency in Hz
end
if(~exist('SNR','var')||isempty(SNR))
    SNR = 10.5; %SNR in dB
end
if(~exist('ncpi','var')||isempty(ncpi))
    ncpi = 1; %SNR in dB
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
if(~exist('compute_bandgraph','var')||isempty(compute_bandgraph))
    compute_bandgraph = false;
end
indexs = 1;
% physical parameters
c = 299792458; % speed of light in m/s
k  = 1.3806503*10^-23; % Boltzmann's constant

% calculate beamsize in degrees
wavelength = c/operating_frequency;
beamfactor=17/20; %factor to convert to 3dB beamwidth

azbw_transmit = asind(wavelength/width) * beamfactor ;
elbw_transmit = asind(wavelength/height) * beamfactor ;

azbw = asind(wavelength/width) * beamfactor;
elbw = asind(wavelength/height) * beamfactor;


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
    elbw_transmit = asind(wavelength/(height * spoil_parameter))  * beamfactor ;


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
            beam_el_transmit(count) = beamcenter_el;

            beam_ranges_transmit(count)=sargmtiangles2(radar_height,beam_el_transmit(count));
            beam_max_range_transmit(count) = sargmtiangles2(radar_height,beam_el_transmit(count) + elbw_transmit/2);
            if( abs( beam_max_range_transmit(count))~= beam_max_range_transmit(count))
                beam_max_range_transmit(count) = max_range;
            end
            if(beam_max_range_transmit(count)>max_range)
                beam_max_range_transmit(count) = max_range;
            end
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
            if(abs(beam_max_range_transmit(count))~= beam_max_range_transmit(count))
                beam_max_range_transmit(count) = max_range;
            end
            if(beam_max_range_transmit(count)>max_range)
                beam_max_range_transmit(count) = max_range;
            end
            if(~isreal(beam_max_range_transmit(count)))
                beam_max_range_transmit(count) = max_range;
            end
            beam_x_transmit(count) = beam_ranges_transmit(count) * sind(beam_az_transmit(count));
            beam_y_transmit(count) =  beam_ranges_transmit(count) * cosd(beam_az_transmit(count));
            if(length(azbw_transmit)<2)
                PSG_transmit = 10*log10(4 * pi * width2*height/wavelength^2);
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
doppler_dwell_time_beam  = doppler_dwell_time_beam  * ncpi;
%fprintf(1,'TOTAL DWELL TIME NEEDED FOR DOPPLER %3.1f\n',sum(total_doppler_dwell_time));




% %% beam power
% for ii = 1:length(beam_max_range)
%     Ptav(ii) = (4 * pi)^3 * (beam_max_range(ii) * 1000).^4 * 10^(NF/10) * k * T * 10^(losses/10) * 10^(SNR/10)./...
%                     (10.^(PSG_transmit/10)  * doppler_dwell_time_beam(ii) *10^(PSG/10).* cosd(beam_az(ii)).^2 .* cosd(beam_el(ii)).^2);
% end
% compute dwell time based on power
beam_az_transmit(find(abs(beam_az_transmit)>60) ) = 60;

power_dwell_time = (4 * pi)^2 * (beam_max_range_transmit * 1000).^4 * 10^(NF/10) * k * T * 10^(losses/10) * 10^(SNR/10)./...
                    (10.^(PSG_transmit/10) * Ptransmit *   rcs * width * height .* cosd(beam_az_transmit).^3 .* cosd(beam_el_transmit).^3);
power_dwell_time = power_dwell_time * ncpi;                
%fprintf(1,'TOTAL DWELL TIME NEEDED FOR POWER %3.1f\n',sum(power_dwell_time));
%5
   power_dwell_time_band = [];
    doppler_dwell_time_band = [];
    X = [];
    Y = [];
    res = 1;
    if(compute_bandgraph)
        count1 = 1;
        count2 =  1;

        
        for range = min_range:res:max_range
            for az = [azimuth_search(1):res:azimuth_search(2)]
                range_all(count2,count1) = range;
                X(count2,count1) = range * sind(az);
                Y(count2,count1) = range * cosd(az);
                [el,grazing]=sargmtiangles(radar_height,range);
                range_from_beam = sqrt( (X(count2,count1) - beam_x_transmit).^2 +  (Y(count2,count1) - beam_y_transmit).^2 );
                [j closest_beam]  = min(range_from_beam);
                closest_beam = closest_beam(1);
                if(length(PSG_transmit)==1)
                    closest_beam = 1;
                end
                power_dwell_time_band(count2,count1) =  (4 * pi)^2 *  (range  * 1000).^4 * ...
                    10^(NF/10) * k * T * 10^(losses/10) * 10^(SNR/10)./...
                    (10.^(PSG_transmit(closest_beam)/10) * Ptransmit *   rcs *...
                    width * height .* cosd(az).^2 .* cosd(el).^2);

                doppler_dwell_time_band(count2,count1)  = wavelength/min_doppler_resolution * 1./sind(90-max(abs(el)));
               % receive_beams_band(count2,count1) = recv_beams(closest_beam);
                count1 = count1 + 1;
            end
            count1 = 1;
            count2 = count2 + 1;
        end
    end
    doppler_dwell_time_band = doppler_dwell_time_band * ncpi;
    power_dwell_time_band = power_dwell_time_band * ncpi;