Fs = 600E3; %samples/second
Fc = 1.3E9; %hertz
C = 299792458; %m/s
prf = 2000; % in hz
pri = 1/prf;
pulsSamples = pri * Fs;
pulseBandwidth = 300E3; % in hz
dutyFactor = .9;
minimumRange = C*pri*dutyFactor;
pulseLength = pri*dutyFactor ; % in seconds
gamma = -10;
%geometry of transmit array
% [geometry w h]= loadgeometry('rectangle',[],arrayLength,arrayHeight,lamda);
% noise figure of transmitter in dB
NF = 4;
% losses in dB
losses = 5;
% receiver temperature in Kelvin
T = 290;
% Boltzmann's constant
k  = 1.3806503*10^-23;
gainTransmit = 0;
gainReceive = 0;
% transmit power in watts
powerTransmit = 100/dutyFactor;
invert_image = true;
lamda = C/Fc ;
% width of array in meters
width = 12*lamda;
% height of array in meters
height = lamda;
% constant gamma cluter
gamma = -10;
% 3:1 overlap (75%)
overlap = 0;
[fullarrayGeometry]= load_geometry('line',[],width,height,299792458/Fc );
fullarrayGeometry = fullarrayGeometry*.5;
x = fullarrayGeometry(:,2)
fullarrayGeometry(:,2) = fullarrayGeometry(:,1);
fullarrayGeometry(:,1) = x;
lamda = C/Fc;
snr = radarequation(k*T,losses,NF,powerTransmit,gainTransmit,gainReceive,Fs,lamda);
% derived parameters
lamda = C/Fc;

% dwell_time = 1/(rpm/60) * horizontalBeamwidth/360;
Npulses = floor(time * prf);%number of pulse in look direction
plength = round(Fs/ prf);