
transmitters = 1;
%center frequency in Hz
operating_frequency = 5E9;
% width of array in meters
width = 1.2;
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
Npulses = 30;
%length of pulse in seconds
pulse_length = 50E-6;

% sampling rate
sampling_rate = 5E6/4;
% bandwidth of chirp
chirp_bandwidth = 4E6/4;
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
Ptransmit = 200;
% Boltzmann's constant
k  = 1.3806503*10^-23;
pre_doppler_response = [0 0 1 0 0 ];
prf = 1500;
make_targets = false;
make_clutter = true;
make_noise = false;
compute_stap = false;
compute_radarequation = false;
compute_subarray = false;
line_array = true;
range_gate = 340;
az = [-10:10];
el = [-10:10];
range_ambiguous = false;
compute_post_doppler_stap = false;
spacing = .12;

% power clipping in dB on receiver
power_clipping = 50;
%% create target and clutter parameters
number_of_targets = 1;
% range of ranges, in m, to gmti search
minmaxrange =  [20 100];
%radar velocity in m/s
radar_velocity = [18; 0; 0];
%target position in m
radar_position= [0; 0; 5];
%range of clutter in m (from start to horizon)
rearth=6446*4/3;
clutter_range = [5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];

% create targets
[target_position target_velocity] = create_targets(number_of_targets);
%% look direction
% transmit and receive azimuth and elevation in degrees
% transmit and receive azimuth and elevation in degrees
transmit_beam = [0 -4.9];
receive_beam  = [0 -4.9];

%dwell time in seconds
dwell_time = .1;

