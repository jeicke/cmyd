fprintf(1,'***Starting simulation\n');
% number of channels in geometry
sims.channels = 10;
% number of targets to simulate
sims.number_of_targets = 1;
% range window to simulate (the maximum and minimum range in simulation)
sims.range_window = [1 100];
% range of clutter 
sims.clutter_range = [25 35];
% zero means no variation over range, anything else is the magnitude of the
% random clutter variation over range
sims.clutter_range_snr = [0];
sims.velocity = [100 0 0];
sims.altitude = 10 * 1000;
% jammer SNR
sims.JNR = 40;

% radar related parameters
radar.range = [5 100];
radar.pulse_length = 100E-6;
radar.Npulses = 100;
radar.prf = 300;
radar.sampling_bandwidth = 5E6;
radar.operating_frequency = 400E6;
radar.chirp_bandwidth = 4E6;
radar.c = 299792458;
fprintf(1,'***Reading Geometry\n');
radar.geometry = load_geometry('line',sims.channels);

% generate chirp
fprintf(1,'***Generating Chirp\n');
radar.pulse_replica = m_chirp(radar.chirp_bandwidth,radar.pulse_length,-radar.chirp_bandwidth/2,radar.sampling_bandwidth);


% create power as a function of angle
resolution = 1;
clutter_power = create_power(radar,resolution,[],[],sims);
%create doppler for each clutter ring
dopplers = create_doppler(radar,resolution,sims);


% create clutter
fprintf(1,'***Creating Clutter\n');
EL = 0.0;
clutter = create_clutter(clutter_power,dopplers,(radar.operating_frequency)/299792458,radar.geometry,resolution,EL,radar.Npulses,radar.prf);
D=clutter;
[D radar] = pulse_compression(D,radar,'cheb',70);
[D radar] = doppler_filter(D,radar,128,'cheb',70);
bearings = [270];
[V] = steering_vector(radar.geometry, (bearings ) * pi/180, ...
    zeros(size(bearings )), (radar.operating_frequency)/299792458 );
[D radar] = beamform(D,radar,1:size(D,3),V);
a = 1;


