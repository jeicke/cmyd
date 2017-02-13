%center frequency in Hz
operating_frequency = 5E9;
% width of array in meters
width = .6;
% height of array in meters
height = .4;
% constant gamma cluter
gamma = -10;
% 3:1 overlap (75%)
overlap = .75;
% speed of light in m/s
c = 299792458;
subarray_sidelobes = 40;
%subarrray elements in each diminsion
subarray_elements_width = 12;
subarray_elements_height = 12;
% pulses in pulse doppler radar
 Npulses = 64;
%length of pulse in seconds
pulse_length = 50E-6;

% sampling rate 
sampling_rate = 5E6/8;
% bandwidth of chirp
chirp_bandwidth = 4E6/8;
% resolution in degrees of clutter model
resolution = .125;
power_threshold = 40;
%
NF = 3;
% losses in dB
losses = 5;
% receiver temperature in Kelvin
T = 290;
% transmitted power in Watts
Ptransmit = 100;
% Boltzmann's constant
k  = 1.3806503*10^-23;
%% create geometry of radar

% receive subarray geometry
[subarray_geometry]= load_geometry('subarray',[],subarray_elements_width,subarray_elements_height,299792458/5E9);
% width and height of subarray in meters
subarray_width =  max(abs(subarray_geometry(:,1))) * 2;
subarray_height =  max(abs(subarray_geometry(:,3))) * 2;

% gain of subarray
PSG_subarray = 10*log10(4 * pi * subarray_width *subarray_height /(299792458/operating_frequency )^2);

% gain of full array
PSG = 10*log10(4 * pi * width * height /(299792458/operating_frequency )^2);
%subarray weighting
w1 = chebwin(subarray_elements_width,subarray_sidelobes );
w2 = chebwin(subarray_elements_height,subarray_sidelobes );
subarray_window = kron(w1/sum(w1),w2/sum(w2));

% now create full receive array geometry
[fullarray_geometry w h]= load_geometry('rectangle',[],width,height,2 * (subarray_width * (1-overlap)));

%transmit array geometry

[transmit_geometry w h]= load_geometry('rectangle',[],width,height,299792458/operating_frequency);
%% create target and clutter parameters
number_of_targets = 1;
% range of ranges, in m, to gmti search
minmaxrange =  [5 100];
%radar velocity in m/s
radar_velocity = [18 * 4; 0; 0];
%target position in m
radar_position= [0; 0; 5];
%range of clutter in m (from start to horizon)
rearth=6446*4/3;
clutter_range = [.1 sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];

% create targets
[target_position target_velocity] = create_targets(number_of_targets);
prf = 6000;%c/(2 * (clutter_range (2)+5) * 1000);%ceil((350) * 4 * operating_frequency/c);%6 * norm(radar_velocity)/2.4;

%% look direction
% transmit and receive azimuth and elevation in degrees
transmit_beam = [0 -2];
receive_beam  = [0 -2];
%dwell time in seconds
dwell_time = .1;
%% derived parameters
    unambiguous_range = c/(2 * prf)/1000;

    pulse_replica = m_chirp(chirp_bandwidth,pulse_length,-chirp_bandwidth/2,sampling_rate );
    % create targets
    
    % range of Azimuths to simulate for clutter
     AZ = single( ([-90:resolution:90]) * pi/180);

   
    % calculate range spacing in meters
    rngspa= 299792458/(sampling_rate)/2;
    
    %ranges from start of recording to horizon
    ranges = [100:rngspa:(clutter_range(2))*1000];
    
    %set of ranges that occur between pulses
    unambiguous_range = [100:rngspa:(unambiguous_range) * 1000];

%% simulation
    fprintf(1,'***Creating Targets\n');
    %create target
    range_ambiguous = false;
    target = target_generate(target_position,radar_position,radar_velocity-target_velocity,...
        pulse_replica,ranges , sampling_rate,transmit_beam ,...
        k*T,losses,dwell_time,NF,1,Ptransmit,PSG ,operating_frequency,Npulses,prf,...
        transmit_geometry,fullarray_geometry,subarray_geometry,subarray_window,PSG_subarray ,unambiguous_range,range_ambiguous);
    
    fprintf(1,'***Creating Clutter\n');
    % clutter power as a function of angle
    clutter_power = create_power(AZ,pulse_replica,ranges ,sampling_rate,transmit_beam,...
        radar_position(3),k*T,losses,dwell_time,NF,c/operating_frequency,gamma,Ptransmit,...
        PSG,operating_frequency,fullarray_geometry,transmit_geometry,subarray_geometry,subarray_window,PSG_subarray );
    %doppler for each clutter ring
    [dopplers EL]= create_doppler(AZ,radar_velocity, ranges,operating_frequency,radar_position(3));

    % create clutter
    
     [clutter]= create_clutter4(AZ,clutter_power,dopplers,(operating_frequency)/299792458,fullarray_geometry,EL,Npulses,prf,ranges,unambiguous_range, power_threshold,range_ambiguous);
    
    % add noise
     fprintf(1,'***Adding noise \n');
    noise = (1/sqrt(2))*(randn(size(target)) + 1i * randn(size(target)));
    % compute receive beam
    D=  noise  * sqrt(pulse_length * chirp_bandwidth) *  sqrt(Npulses) + target + clutter;
    D = clutter;
    %% process radar
     fprintf(1,'***Processing radar\n');
    [Ds ] = pulse_compression(D,pulse_replica,'cheb',70);
    Ds = Ds(length(pulse_replica):end,:,:);
    [D] = doppler_filter(Ds,Npulses,'none',70,'false');


    [V] = steering_vector(fullarray_geometry, receive_beam(1) * pi/180, ...
         receive_beam(2) * pi/180, (operating_frequency)/299792458 );
    [Dbf ] = beamform(D,1:size(D,3),V);
    
     az = [-20:20];
    f_simo = prf * (-size(Dbf,2)/2:size(Dbf,2)/2-1)/(size(Dbf,2)/2);
    [V] = steering_vector(fullarray_geometry, az * pi/180, ...
         receive_beam(2) * pi/180 * ones(size(az)), (operating_frequency)/299792458 );
    [Dbf_simo ] = beamform(D,1:size(D,3),V );
  
    %     p1 = 20*log10(abs(D(19,:,1)));
    %     plot(p1 -median(p1(1:30)),'r')
    %     hold on;
    %     p2 = 20*log10(abs(Dbf(19,:,1)));
    %     plot(p2 -median(p2(1:30)),'c')
    
    % plot
    rngspa= 299792458/(sampling_rate)/2;


   % fprintf(1,'Steering Transmit %3.1fX%3.1f Steering Receive %3.1fX%3.1f \n',beam_transmit(1,beam),beam_transmit(2,beam),beam_receive(1,beam),beam_receive(2,beam))
   
    p = 20*log10(abs(Dbf(:,:,1)));
    unambiguous_range = c/(2 * prf)/1000;
    unambiguous_range = [length(pulse_replica) * rngspa:rngspa:(unambiguous_range) * 1000];
    imagesc(p)
    caxis([0 90])
    colorbar
    drawnow;
return;
    fprintf(1,'***Processing radar STAP\n');
     [V] = steering_vector(fullarray_geometry, receive_beam(1) * pi/180, ...
         EL(400), (operating_frequency)/299792458 );
     
    [Dbf2 w] = beamform_stap(Ds,V'.',[0 1 0]',4,50,1E-4,4);
    
    %[sinr sinr2] = calculate_sinr(w,  clutter,V,Npulses,size(fullarray_geometry,1),prf ,20,20,20);
    
    [Dbf2 shade_doppler] = doppler_filter(Dbf2,Npulses,'cheb',70);
    [Dbf2] = mti_normalize(Dbf2,'median');
    p = 20*log10(abs(Dbf2(:,:,1)));
    figure;
    imagesc(v, unambiguous_range/1000,p)
    caxis([0 70])
    colorbar
    drawnow;