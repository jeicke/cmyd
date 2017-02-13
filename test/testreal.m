% Test of proper scattering range
% tests that we have the poper range to scatterer
% RADAR PARAMETERS

Fs = 1024; %samples/second
Fc = 0; %hertz
rpm = 24; %rotations per minute of radar
verticalBeamwidth = 25; % in degrees
horizontalBeamwidth = 1.8; % in degrees
prf = 1000; % in hz
pulseBandwidth = 0; % in hz
pulseLength =  0.0977; % in seconds
C = 600; %m/s
filebase = 'C:\Users\rnl\Desktop\testdata\testmulti';
image_filename = 'C:\Users\rnl\Desktop\ISAR code\test3.png';
invert_image = true;
mask;
% source and receiver locations
receiverPosition = [0 -200 0];   
transmitterPosition = [100 0 0];
receiverVelocity = 10 * [0 1 0];   
transmitterVelocity= [0 0 0];

% derived parameters
lamda = C/Fc;
arrayHeight = (1.27 * lamda)/(verticalBeamwidth *pi/180);
arrayLength = (1.27 * lamda)/(horizontalBeamwidth *pi/180);
dwell_time = 1/(rpm/60) * horizontalBeamwidth/360;
Npulses = floor(dwell_time * prf);%number of pulse in look direction
plength = round(Fs/ prf);
% just do a few look for this simulation
time = 40;

% set track for rotating transmit antenna
t = [0:dwell_time:time];
orientation = zeros(4,length(t));
orientation(3,:) = [0:horizontalBeamwidth:(length(t)-1)*horizontalBeamwidth],360;% - radar_velocity(3)/2  * time;
orientation(4,:) = t;

%geometry of transmit array
[geometry w h]= loadgeometry('rectangle',[],arrayLength,arrayHeight,lamda);

% make transmitter and receiver tracks
[ track] =maketrack(time,.1,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );

% make scattered point pegged to transmitter (like we're on a boat)



clear reflextivity;
for ii = 1:size(scatteredpoint,2);
    %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    
    reflextivity(ii) = sigma(ii);
    track(:,:,ii+2) = track(:,:,1);
    track(1,:,ii+2) =track(1,:,ii+2) +scatteredpoint(1,ii);
    track(2,:,ii+2) =track(2,:,ii+2) +scatteredpoint(2,ii);
    track(3,:,ii+2) =track(3,:,ii+2) +scatteredpoint(3,ii);
end
w =chebwin(40,70);
w = kron(w,[1 1 1]');
w = w/sum(abs(w));
geometry = [];
w = [];
radsim(filebase,time,track,0,...
       'timeseries' ,@chirppulse,'parameters',[30 pulseLength prf],...
       'fs',Fs,'fc',Fc,'transmitgeometry',geometry,...
       'transmitorientation',orientation,'clutter',reflextivity',...
       'transmitshading',w,'c',C); 