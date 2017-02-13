%% parameters for radar
%center frequency in Hz
operating_frequency = 5E9;
% transmitted power in Watts
Ptransmit = 100;
% width of array in meters
width = 2.4;
% height of array in meters
height = .4;
%receiver noise figure in dB
NF = 3;
% losses in dB
losses = 5;
% receiver temperature in Kelvin
T = 290;
% azimuthal scan window in start and end degrees
azimuth_search = [0 0];
% transmit spoil in fractional amount
spoil_parameter = 1;
% min detectable velocity in m/s
min_detectable_velocity = 1;
% receive beams
sim_recv_beams = 1;
% set for explicit beam ranges
elev_beam_ranges = [];
% amount to spoil beam
spoil = [];
% radar update rate
radar_update_rate = 1;
% constant gamma cluter
gamma = -10;
% 3:1 overlap (75%)
overlap = .75;
% speed of light in m/s
c = 299792458;
% Boltzmann's constant
k  = 1.3806503*10^-23;
%sub array sidelobes
subarray_sidelobes = 40;
%subarrray elements in each diminsion
subarray_elements_width = 12;
subarray_elements_height = 12;
% pulses in pulse doppler radar
Npulses = 32;
%length of pulse in seconds
pulse_length = 50E-6;

% sampling rate 
sampling_rate = 5E6/4;
% bandwidth of chirp
chirp_bandwidth = 4E6/4;
% resolution in degrees of clutter model
resolution = .125;
power_threshold = 40;
%% create geometry of radar

% receive subarray geometry
[subarray_geometry]= load_geometry('subarray',[],subarray_elements_width,subarray_elements_height,299792458/5E9);
% width and height of subarray in meters
subarray_width =  max(abs(subarray_geometry(:,1))) * 2;
subarray_height =  max(abs(subarray_geometry(:,3))) * 2;

% gain of subarray
PSG_subarray = 10*log10(4 * pi * subarray_width *subarray_height /(299792458/operating_frequency )^2);

%subarray weighting
w1 = chebwin(subarray_elements_width,subarray_sidelobes );
w2 = chebwin(subarray_elements_height,subarray_sidelobes );
subarray_window = kron(w1/sum(w1),w2/sum(w2));

% now create full receive array geometry
[fullarray_geometry w h]= load_geometry('rectangle',[],width,height,2 * (subarray_width * (1-overlap)));

%transmit array geometry

[transmit_geometry w h]= load_geometry('rectangle',[],width,height,299792458/operating_frequency);
%
% HGA.theta_source = 0;
% HGA.phi_source = 0;
% HGA.SLx = fullarray_geometry;
% HGA.subarray_geometry = subarray_geometry ;
% HGA.subarray_theta_source =  0;
% HGA.subarray_phi_source = 0;
% HGA.subarray_window = subarray_window ;
% [theory_matrix az de] = calcBeamPattern(5E9,[-90:.1:90],0,[],HGA);
% plot(az,10*log10(theory_matrix),'g')
%% create target and clutter parameters
number_of_targets = 1;
% range of ranges, in m, to gmti search
minmaxrange =  [20 100];
%radar velocity in m/s
radar_velocity = [0; 125; 0];
%target position in m
radar_position= [0; 0; 5];
%range of clutter in m (from start to horizon)
rearth=6446*4/3;
clutter_range = [.1 sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];

% create targets
[target_position target_velocity] = create_targets(number_of_targets);

%% create grid of beams for target based on radar parameters
[beam_transmit beam_receive doppler_dwell_time_beam ] = gmtigrid(operating_frequency,[],Ptransmit,...
    radar_position(3),minmaxrange,width,height,NF,losses,...
    T,[],azimuth_search,spoil_parameter,min_detectable_velocity,...
    sim_recv_beams,elev_beam_ranges,spoil);


%% step through transmit beams

for beam = 1:size(beam_transmit,2)
    fprintf(1,'***Working on transmit beam %d of %d\n',beam,size(beam_transmit,2));
    beam = 1;
    % creat beam parameters
     %beam_transmit(:,beam) = [90; 0;1;minmaxrange(2);minmaxrange(1);35;.01];
     %beam_receive(:,beam) = [90; 0;1;minmaxrange(2);minmaxrange(1);35];

    % 3 times main beam clutter width
    %c/(2 * (minmaxrange(2)+5) * 1000);
    prf = 3000;%c/(2 * (clutter_range (2)+5) * 1000);%ceil((350) * 4 * operating_frequency/c);%6 * norm(radar_velocity)/2.4;
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
    
    %ranges covered in receive beam
    max_range =  [(beam_transmit(5,beam))*1000:rngspa:(beam_transmit(4,beam))*1000];
    fprintf(1,'***Creating Targets\n');
    %create target
    target = target_generate(target_position,radar_position,radar_velocity-target_velocity,...
        pulse_replica,ranges , sampling_rate,[beam_transmit(1,beam) beam_transmit(2,beam)],...
        k*T,losses,beam_transmit(7,beam),NF,1,Ptransmit,beam_transmit(6,beam),operating_frequency,Npulses,prf,...
        transmit_geometry,fullarray_geometry,subarray_geometry,subarray_window,PSG_subarray ,unambiguous_range,max_range);
    
    fprintf(1,'***Creating Clutter\n');
    % clutter power as a function of angle
    clutter_power = create_power(AZ,pulse_replica,ranges ,sampling_rate,[beam_transmit(1,beam) beam_transmit(2,beam)],...
        radar_position(3),k*T,losses,beam_transmit(7,beam),NF,c/operating_frequency,gamma,Ptransmit,...
        beam_transmit(6,beam),operating_frequency,fullarray_geometry,transmit_geometry,subarray_geometry,subarray_window,PSG_subarray );
    %doppler for each clutter ring
    [dopplers EL]= create_doppler(AZ,radar_velocity, ranges,operating_frequency,radar_position(3));

    % create clutter
    
     %clutter_power = clutter_power(:,500);
    % dopplers = dopplers(:,500);
%     
     [clutter SVa]= create_clutter(AZ,clutter_power,dopplers,(operating_frequency)/299792458,fullarray_geometry,EL,Npulses,prf,ranges,unambiguous_range, power_threshold);
    % trim clutter patches
%     x = 10*log10(max(abs(clutter_power),[],2));
%     [c ii] = find(x>(max(max(x))-power_threshold));
%     clutter_power = clutter_power(c,:);
%     dopplers = dopplers(c,:);
%     AZ = AZ(c);
   
    % add noise
     fprintf(1,'***Adding noise \n');
    noise = (1/sqrt(2))*(randn(size(target)) + 1i * randn(size(target)));
    % compute receive beam
    D=  noise  * sqrt(pulse_length * chirp_bandwidth) *  sqrt(Npulses) + target + clutter;
    %D = clutter;

    % process radar
     fprintf(1,'***Processing radar\n');
    [Ds ] = pulse_compression(D,pulse_replica,'cheb',70);
    Ds = Ds(length(pulse_replica):end,:,:);
    [D] = doppler_filter(Ds,64,'cheb',70);


    [V] = steering_vector(fullarray_geometry,  beam_receive(1,beam) *pi/180, ...
        beam_receive(2,beam)* pi/180, (operating_frequency)/299792458 );
    [Dbf ] = beamform(D,1:size(D,3),V'.' );
    %     p1 = 20*log10(abs(D(19,:,1)));
    %     plot(p1 -median(p1(1:30)),'r')
    %     hold on;
    %     p2 = 20*log10(abs(Dbf(19,:,1)));
    %     plot(p2 -median(p2(1:30)),'c')
    
    % plot
    rngspa= 299792458/(sampling_rate)/2;
    ranges = [beam_transmit(5,beam)*1000:rngspa:beam_transmit(4,beam)*1000];

    fprintf(1,'Steering Transmit %3.1fX%3.1f Steering Receive %3.1fX%3.1f \n',beam_transmit(1,beam),beam_transmit(2,beam),beam_receive(1,beam),beam_receive(2,beam))
    v = prf * c/operating_frequency/4 * ([-Npulses/2:Npulses/2]/Npulses/2);
    p = 20*log10(abs(Dbf(:,:,1)));
    unambiguous_range = c/(2 * prf)/1000;
    unambiguous_range = [length(pulse_replica) * rngspa:rngspa:(unambiguous_range) * 1000];
    imagesc(p)
    caxis([0 90])
    colorbar
    drawnow;

    %    Dbf2 = Dbf;
     %   bearings2 = [-180:2:180];
     fprintf(1,'***Processing radar STAP\n');
    [V] = steering_vector(fullarray_geometry, beam_receive(1,beam) *pi/180, ...
        EL(1100), (operating_frequency)/299792458 );
    [Dbf2 w] = beamform_stap(Ds,SVa'.',[0 1 0]',20,10,1E-4,20);
    
    sinr = calculate_sinr(w,  clutter,SVa'.',Npulses,size(fullarray_geometry,1),prf ,10,20,20);
    %sinr = calculate_sinr(w,clutter,V'.',Npulses,size(fullarray_geometry,1),20,20,20);
   % T=toeplitz([1 -2 1 zeros(1,pulses-3)]',[1 zeros(1,size(D,2)-1)]);
    [Dbf2 shade_doppler] = doppler_filter(Dbf2,100,'cheb',70);
    [Dbf2] = mti_normalize(Dbf2,'median');
    p = 20*log10(abs(Dbf2(:,:,1)));
    figure;
    imagesc(v, unambiguous_range/1000,p)
    caxis([0 90])
    colorbar
    drawnow;
    break;
end

