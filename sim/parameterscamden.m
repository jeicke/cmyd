%make bistatic
%
ISBISTATIC = true;
ISBISTATIC = false;
transmitters = 1;
c = 299792458;
%center frequency in Hz
operating_frequency = 1.3E9;
lamda = c/operating_frequency ;
% width of array in meters
width = 6*lamda;
% height of array in meters
height = lamda;
% constant gamma cluter
gamma = -10;
% 3:1 overlap (75%)
overlap = 0;
% speed of light in m/s

subarray_sidelobes = 40;
%subarrray elements in each diminsion
subarray_elements_width = 1;
subarray_elements_height = 1;
% pulses in pulse doppler radar
Npulses = 10;
prf = 1E3;
%length of pulse in seconds
pulse_length = 1/prf * .1;

% sampling rate
sampling_rate = 200E3;
% bandwidth of chirp
chirp_bandwidth = 100E3;
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
prf = 2000;
make_targets = false;
make_clutter = true;
make_noise = false;
compute_stap = false;
compute_radarequation = false;
compute_subarray = false;
line_array = true;
range_gate = 1;
az = [-180:180];
el = [-180:180];
range_ambiguous = false;
compute_post_doppler_stap = false;


% power clipping in dB on receiver
power_clipping = 50;
%% create target and clutter parameters
number_of_targets = 1;
% range of ranges, in km, to gmti search
minmaxrange =  [12.5 25];
%radar velocity in m/s
radar_velocity = [60; 0; 0];
%target position in km
radar_position= [0; 0; 10];
% if is bistatic
if(ISBISTATIC)
    %target position in km
    emitter_position= [-32.5; 0; 20];
end
%range of clutter in m (from start to horizon)
rearth=6446*4/3;
clutter_range = [12.5 25];%[5  sqrt(2 * rearth * radar_position(3) + radar_position(3)^2)];

% create targets
[target_position target_velocity] = create_targets(number_of_targets);
%% look direction
% transmit and receive azimuth and elevation in degrees
% transmit and receive azimuth and elevation in degrees
transmit_beam = [0 -4.9];
receive_beam  = [0 -4.9];

%dwell time in seconds
dwell_time = .1;

