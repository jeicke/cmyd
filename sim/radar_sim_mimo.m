%center frequency in Hz
operating_frequency = 5E9;
% width of array in meters
width = 2.4;
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
Npulses = 128;
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

transmitters = 2;
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

clear mimo_transmit_geometry;
mimo_transmit_geometry(:,:,1) = [2* mean(fullarray_geometry(1: size(fullarray_geometry,1)/2,1)) 0 0 ];
mimo_transmit_geometry(:,:,2) = [2* mean(fullarray_geometry( size(fullarray_geometry,1)/2+1: size(fullarray_geometry,1),1)) 0 0 ];

mimo_receive_geometry = [fullarray_geometry; fullarray_geometry];
mimo_receive_geometry(1:size(fullarray_geometry,1)/2,1) = mimo_receive_geometry(1:size(fullarray_geometry,1)/2,1) + 2*mean(fullarray_geometry(1:size(fullarray_geometry,1)/2,1));
mimo_receive_geometry(size(fullarray_geometry,1)/2+1:size(fullarray_geometry,1),1) = mimo_receive_geometry(size(fullarray_geometry,1)/2+1:size(fullarray_geometry,1),1) + 2*mean(fullarray_geometry( size(fullarray_geometry,1)/2+1: size(fullarray_geometry,1),1));
%% create target and clutter parameters
number_of_targets = 1;
% range of ranges, in m, to gmti search
minmaxrange =  [5 100];
%radar velocity in m/s
radar_velocity = [72; 0; 0];
%target position in m
radar_position= [0; 0; 5];
%range of clutter in m (from start to horizon)
rearth=6446*4/3;
clutter_range = [.1 sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];

% create targets
[target_position target_velocity] = create_targets(number_of_targets);

prf = 12000;%c/(2 * (clutter_range (2)+5) * 1000);%ceil((350) * 4 * operating_frequency/c);%6 * norm(radar_velocity)/2.4;

%% look direction
% transmit and receive azimuth and elevation in degrees
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
target = target_generate_mimo2(target_position,radar_position,radar_velocity-target_velocity,...
    pulse_replica,ranges , sampling_rate,transmit_beam ,...
    k*T,losses,dwell_time,NF,1,Ptransmit,PSG ,operating_frequency,Npulses,prf,...
    transmit_geometry,fullarray_geometry,subarray_geometry,subarray_window,mimo_transmit_geometry,PSG_subarray ,unambiguous_range,transmitters,range_ambiguous);

fprintf(1,'***Creating Clutter\n');
% clutter power as a function of angle
   clutter_power = create_power(AZ,pulse_replica,ranges ,sampling_rate,transmit_beam,...
        radar_position(3),k*T,losses,dwell_time,NF,c/operating_frequency,gamma,Ptransmit,...
        PSG,operating_frequency,fullarray_geometry,transmit_geometry,subarray_geometry,subarray_window,PSG_subarray );
%doppler for each clutter ring
[dopplers EL]= create_doppler(AZ,radar_velocity, ranges,operating_frequency,radar_position(3));

% create clutter

%clutter_power = clutter_power(:,500);
% dopplers = dopplers(:,500);
%

[clutter]= create_clutter_mimo(AZ,clutter_power,dopplers,(operating_frequency)/299792458,fullarray_geometry,mimo_transmit_geometry,EL,Npulses,prf,ranges,unambiguous_range, power_threshold, transmitters,transmit_beam ,range_ambiguous );
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
%D=  noise  * sqrt(pulse_length * chirp_bandwidth) *  sqrt(Npulses) + target + clutter;
D =  clutter;
%% process radar
% process radar
fprintf(1,'***Processing radar\n');
[Ds ] = pulse_compression(D,pulse_replica,'cheb',70);
Ds = Ds(length(pulse_replica):end,:,:);
[D] = doppler_filter(Ds,Npulses,'none',70,false);
clear D4;
D4(:,:,1:1:size(fullarray_geometry,1)) = D(:,1:(Npulses/2),:);
D4(:,:,size(fullarray_geometry,1)+1:2*size(fullarray_geometry,1)) = D(:,(Npulses/2+1):Npulses,:);
%[D4] = inverse_doppler_filter(D4);
%[D4] = doppler_filter(D4,Npulses,'cheb',70);


fprintf(1,'***Processing radar STAP\n');
 az = [0];
 [V1] = steering_vector(mimo_receive_geometry, az * pi/180, ...
     receive_beam(2) * pi/180 * ones(size(az)), (operating_frequency)/299792458 );
% 
[Dbf_mimo ] = beamform(D4,1:size(D4,3),V1 );
Npulses = Npulses/2;
f_mimo = prf/2 * (-size(Dbf_mimo,2)/2:size(Dbf_mimo,2)/2-1)/(size(Dbf_mimo,2)/2);
imagesc(fftshift((20*log10(abs(Dbf_mimo(:,:,1)))),2))
caxis([0 90])
    colorbar
    drawnow;
clear clutter;
clear target;
clear noise;
return;
% rngspa= 299792458/(sampling_rate)/2;
% v = prf  * c/operating_frequency/4 * ([-Npulses/2:Npulses/2]/Npulses/2);
%  v = prf  * c/operating_frequency/4 * ([0:Npulses-1]/Npulses);
% p = 20*log10(abs(Dbf_mimo(:,:,1)));
% 
% unambiguous_range = c/(2 * prf )/1000;
% unambiguous_range = [length(pulse_replica) * rngspa:rngspa:(unambiguous_range) * 1000];
% imagesc(v,unambiguous_range/1000,p)
% caxis([max(max(p))-70 max(max(p))])
% colorbar
% drawnow;
% [D] = doppler_filter(Ds,Npulses,'cheb',70,false);
% imagesc(10*log10(abs(D(:,:,1))))
% 
% clear D4;
% D4(:,:,1:30) = D(:,1:(Npulses/2),:);
% D4(:,:,31:60) = D(:,(Npulses/2+1):Npulses,:);
% 
% 
% [D] = inverse_doppler_filter(D4);
% 
% % plot(az,10*log10(abs(squeeze(mean(Dbf(400,3,:),1)))),'r',az,10*log10(abs(squeeze(mean(Dbf2(400,3,:),1)))),'g')
% az2 = [0];
% [V] = steering_vector(mimo_receive_geometry, az2 * pi/180, ...
%     EL(400) * ones(size(az2)), (operating_frequency)/299792458 );
% [Dbfs w] = beamform_stap(D,V'.',[0 1 0]',4,5,1E-4,4);
% 
% %[sinr sinr2] = calculate_sinr(w,D,V,Npulses/2,size(mimo_receive_geometry,1),prf/2 ,10,20,20);
% [Dbfs shade_doppler] = doppler_filter(Dbfs,Npulses/2,'none',70);
% [Dbfs] = mti_normalize(Dbfs,'median');
% p = 20*log10(abs(Dbfs(:,:,1)));
% figure;
% imagesc(v, unambiguous_range/1000,p)
% caxis([max(max(p))-70 max(max(p))])
% colorbar
% drawnow;
% 
% 
% % [D] = doppler_filter(Ds,Npulses,'cheb',[],false);
% % D4(:,:,1:30) = D(:,1:(Npulses/2),:);
% % D4(:,:,31:60) = D(:,(Npulses/2+1):Npulses,:);
% %[w v] = beamform_element_stap(D4,V'.',Npulses/2,20,20,20);