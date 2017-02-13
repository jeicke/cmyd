% Test of proper scattering range
% tests that we have the poper range to scatterer
% RADAR PARAMETERS
NUMBEROFSCATTERPOINTS = 50;
Fs = 60E6; %samples/second
Fc = 9410E6; %hertz
rpm = 24; %rotations per minute of radar
verticalBeamwidth = 25; % in degrees
horizontalBeamwidth = 1.8; % in degrees
prf = 3000; % in hz
pulseBandwidth = 30E6; % in hz
pulseLength = 1200E-9; % in seconds
C = 299792458; %m/s

% source and receiver locations
receiverPosition = [0 0 0];   
transmitterPosition = [0 1000 0];
receiverVelocity = [0 0 0];   
transmitterVelocity= [0 0 0];
% derived parameters
lamda = C/Fc;
arrayHeight = (1.27 * lamda)/(verticalBeamwidth *pi/180);
arrayLength = (1.27 * lamda)/(horizontalBeamwidth *pi/180);
dwell_time = 1/(rpm/60) * horizontalBeamwidth/360;
Npulses = floor(dwell_time * prf);%number of pulse in look direction
plength = round(Fs/ prf);
% just do a few look for this simulation
time = dwell_time  * 8;

% set track for rotating transmit antenna
t = [0:dwell_time:time];
orientation = zeros(4,length(t));
orientation(3,:) = [0:horizontalBeamwidth:(length(t)-1)*horizontalBeamwidth];% - radar_velocity(3)/2  * time;
orientation(4,:) = t;

%geometry of transmit array
[geometry w h]= loadgeometry('rectangle',[],arrayLength,arrayHeight,lamda);

% make transmitter and receiver tracks
[ track] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );

% make scattered point pegged to transmitter (like we're on a boat)

clear reflextivity;
radius = 600;
for ii = 1:NUMBEROFSCATTERPOINTS
    scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    
    reflextivity(ii) = 1;
    track(:,:,ii+2) = track(:,:,1);
    track(1,:,ii+2) =track(1,:,ii+2) +scatteredpoint(1);
    track(2,:,ii+2) =track(2,:,ii+2) +scatteredpoint(2);
    track(3,:,ii+2) =track(3,:,ii+2) +scatteredpoint(3);
end
geometry= [];
radsim('C:\Users\rnl\Desktop\testdata\test',time,track,0,...
       'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
       'fs',Fs,'fc',Fc,'transmitgeometry',geometry,'transmitorientation',orientation,'clutter',reflextivity');  
% make replica pulse
pulse_replica=m_chirp(pulseBandwidth,pulseLength,0,Fs)'.';
pulse_replica=pulse_replica(end:-1:1);

% now read data
[data Fs Fc channels] = testread('C:\Users\rnl\Desktop\testdata\test');

% put data into data cube
datapad = [data zeros(1,ceil((length(data)/ plength)/Npulses) *  plength*Npulses-length(data))];
datapad = 100*datapad + 1/sqrt(2) * (randn(size(datapad)) + 1i * randn(size(datapad)));
D = reshape(datapad,plength,Npulses,[]);


% pulse compress data
[Ds ] = pulsecompression(D,pulse_replica,'cheb',70);

% plot track
plottrack(track)
figure;
%plot pulses
rngspa= 299792458/Fs;
ranges = [0:1:size(D,1)-1] * rngspa;
ranges = ranges - (length(pulse_replica)-2) * rngspa;
%%plot result
plot(ranges,(abs(Ds(:,1,1))).^2)
grid on;
xlabel('RANGE (METERS')
ylabel('PULSE COMPRESSED POWER')
trrange = norm(track(1:3,1,1)-track(1:3,1,2));

hold on;
plot([trrange trrange],[0 max((abs(Ds(:,1,1))).^2)],'r--','linewidth',2)
legend('SIMULATION','ACTUAL RANGE')
for ii = 1:(size(track,3)-2)
 trrange = norm(track(1:3,1,1)-track(1:3,1,ii+2));
 trrange = trrange + norm(track(1:3,1,2)-track(1:3,1,ii+2));

 plot([trrange trrange],[0 max((abs(Ds(:,1,1))).^2)],'r--','linewidth',.5)
end