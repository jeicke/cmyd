% Test of proper scattering range
% tests that we have the poper range to scatterer
% RADAR PARAMETERS
NUMBEROFSCATTERPOINTS = 5;
Fs = 60E6; %samples/second
Fc = 9410E6; %hertz
rpm = 24; %rotations per minute of radar
verticalBeamwidth = 25; % in degrees
horizontalBeamwidth = 2.6; % in degrees
prf = 3000; % in hz
pulseBandwidth = 30E6; % in hz
pulseLength = 1200E-9; % in seconds
C = 299792458; %m/s

% source and receiver locations
receiverPosition = [0 0 0];   
transmitterPosition = [0 100 0];
receiverVelocity = [0 0 0];   
transmitterVelocity= [0 100 0];
sourceVelocity= [0 0 0];
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
radius = 800;
for ii = 1:NUMBEROFSCATTERPOINTS
    scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    %scatteredpoint = [100 -200 0]
    reflextivity(ii) = 1;
    track(:,:,ii+2) = track(:,:,1);
    track(1:3,:,ii+2) = 0;
    track(1,:,ii+2) =track(1,:,ii+2) +scatteredpoint(1);
    track(2,:,ii+2) =track(2,:,ii+2) +scatteredpoint(2);
    track(3,:,ii+2) =track(3,:,ii+2) +scatteredpoint(3);
end
geometry= [];
radsim('.\test',time,track,0,...
       'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
       'fs',Fs,'fc',Fc,'transmitgeometry',geometry,'transmitorientation',orientation,'clutter',reflextivity');  
% make replica pulse
pulse_replica=m_chirp(pulseBandwidth,pulseLength,0,Fs)'.';
pulse_replica=pulse_replica(end:-1:1);

% now read data
[data Fs Fc channels] = testread('.\test');

% put data into data cube
datapad = [data zeros(1,ceil((length(data)/ plength)/Npulses) *  plength*Npulses-length(data))];
datapad = 100*datapad + 1/sqrt(2) * (randn(size(datapad)) + 1i * randn(size(datapad)));
D = reshape(datapad,plength,Npulses,[]);


% pulse compress data
[Ds ] = pulsecompression(D,pulse_replica,'cheb',70);
[Df] = dopplerfilter(Ds,128,'cheb',70,false);

% plot track
%plottrack(track,1)
figure;
%plot pulses
rngspa= 299792458/Fs;
ranges = [0:1:size(D,1)-1] * rngspa;
ranges = ranges - (length(pulse_replica)-2) * rngspa;
f = (prf/2)* (-size(Df,2)/2:size(Df,2)/2-1)/(size(Df,2)/2);

%%plot result
%%
figure;
imagesc(f,ranges/1000,20*log10(fftshift((abs(Df(:,:,1))),2)))
hold on;
for ii = 1:(size(track,3)-2)
 p1 = track(1:3,1,1)-track(1:3,1,ii+2);
 p2 = track(1:3,1,ii+2)-track(1:3,1,2);
 t1 = norm(p1);
 t2 = norm(p2);
 u1 = p1/t1;
 u2 = p2/t2;
 
 trrange = t1+t2;
 %p = [x2;y(ii,1);z(ii,1)] /norm( [x2;y(ii,1);z(ii,1)]);
 d1=  (sourceVelocity-transmitterVelocity) * u1* (Fc)/299792458;
 d2=  (receiverVelocity-sourceVelocity) * u2* ( Fc)/299792458;
 doppler = d1 + d2;
 plot([doppler doppler],[ranges(1) ranges(end)]/1000,'r--');
end



p = (track(1:3,1,1)-track(1:3,1,2))/norm((track(1:3,1,1)-track(1:3,1,2)))
doppler = (receiverVelocity-transmitterVelocity)* p * Fc/299792458


plot([doppler doppler],[ranges(1) ranges(end)]/1000,'r')
xlabel('DOPPLER FREQUENCY')
ylabel('RANGE (KM)')