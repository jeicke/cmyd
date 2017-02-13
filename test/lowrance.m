% Test of proper scattering range
% tests that we have the poper range to scatterer
% RADAR PARAMETERS
addpath('C:\Users\jeicke\Documents\programs\msmsmr\ISAR code')
addpath('C:\Users\jeicke\Documents\programs\msmsmr\ISAR code\test')
addpath('C:\Users\jeicke\Documents\programs\msmsmr\ISAR code\bitmaps')
Fs = 300000000; %samples/second
Fc = 9350E6; %hertz
rpm = 24; %rotations per minute of radar
verticalBeamwidth = 25; % in degrees
horizontalBeamwidth = 1.3; % in degrees
% 200-540
prf = 200; % in hz
pulseBandwidth = 75E6; % in hz
pulseLength = 1E-3; % in seconds
C = 299792458; %m/s
filebase = 'D:\testdata\test';
image_filename = 'C:\Users\jeicke\Documents\programs\msmsmr\ISAR code\bitmaps\test5.png';
invert_image = true;
boat_outline;
gamma = -10;
% noise figure of transmitter in dB
NF = 6;
% losses in dB
losses = 5;
% receiver temperature in Kelvin
T = 290;
% Boltzmann's constant
k  = 1.3806503*10^-23;
gainTransmit = 28;
gainReceive = 20;
% transmit power in watts
powerTransmit = 165E-3;

lamda = C/Fc;
snr = radarequation(k*T,losses,NF,powerTransmit,gainTransmit,gainReceive,Fs,lamda);

time =1.2/(rpm/60);
% source and receiver locations
receiverVelocity = 0* [0 1 0];   
receiverPosition = -receiverVelocity * time/2; 
receiverPosition(3) = 2;
transmitterPosition = [-4000 0 16];

transmitterVelocity=[0 0 0];%receiverVelocity;%[0 0 0];

% derived parameters
lamda = C/Fc;
arrayHeight = (1.27 * lamda)/(verticalBeamwidth *pi/180);
arrayLength = (1.27 * lamda)/(horizontalBeamwidth *pi/180);
dwell_time = 1/(rpm/60) * horizontalBeamwidth/360;
Npulses = floor(dwell_time * prf);%number of pulse in look direction
plength = round(Fs/ prf);
% just do a few look for this simulation


% set track for rotating transmit antenna
t = [0:dwell_time:time];
orientation = zeros(4,length(t));
orientation(3,:) = [0:horizontalBeamwidth:(length(t)-1)*horizontalBeamwidth],360;% - radar_velocity(3)/2  * time;
orientation(4,:) = t;

%geometry of transmit array
[geometry w h]= loadgeometry('rectangle',[],arrayLength,arrayHeight,lamda);

% make transmitter and receiver tracks
[ track] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );

% make scattered point pegged to transmitter (like we're on a boat)

%scatteredpoint = scatteredpoint/9;
%scatteredpoint(3,:) = scatteredpoint(3,:) * 2;
clear reflextivity;
for ii = 1:size(scatteredpoint,2);
    %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    
    reflextivity(ii) = sigma(ii);
    track(:,:,ii+2) = track(:,:,1);

    %track(1:3,:,ii+2) = 0;
    track(1,:,ii+2) =track(1,:,ii+2) +scatteredpoint(1,ii);
    track(2,:,ii+2) =track(2,:,ii+2) +scatteredpoint(2,ii);
    track(3,:,ii+2) =track(3,:,ii+2) +scatteredpoint(3,ii);
    
end
snr = 10*log10(snr) *ones(1,length(sigma)+1);
w =hannwindow(length(geometry)/3).';
w = kron(w,[1 1 1]');
w = w/sum(abs(w));
%geometry = [];
 %w = [];
%%
 Parameters = radsim(filebase,time,track,snr,...
       'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
       'fs',Fs,'fc',Fc,'transmitgeometry',geometry,...
       'transmitorientation',orientation,'clutter',reflextivity',...
       'transmitshading',w,'supressdirect',false,'rangecorrect',true ); 
%%
pulse_replica=m_chirp(pulseBandwidth,pulseLength,0,Fs)'.';
 
pulse_replica=pulse_replica(end:-1:1);
rngspa= 299792458/Fs;
[maxDelay minDelay ] = finddelays( track,geometry,C );
minDelay  =  minDelay - 100/Fs;
maxDelay = maxDelay+ 100/Fs;
plength = round(Fs/ prf);
blockLength =plength;
totalLength = blockLength + 2*length(pulse_replica);
ranges = ([1:1:totalLength]-1) * 1/Fs;
ranges =  ranges*C-2*length(pulse_replica)*C/Fs+C/Fs*2;
[h minRangeIndex ] = min(abs(ranges-minDelay *C));
% minRangeIndex  = 3000;
[h maxRangeIndex ] = min(abs(ranges-maxDelay *C));
 
rangeGates = (max([1  minRangeIndex]):min([maxRangeIndex length(ranges)]));
%%
[D blockTimes] = processradar( filebase,prf,1,pulse_replica,rangeGates);

ranges = ranges(rangeGates );

%%
p = squeeze(sum(D,2));
dt = 1/prf;
t = ([0:1:(size(D,3)-1)])*(dt * size(D,2))+dt/2;
azx = interp1(orientation(4,:),orientation(3,:),t);



frameLength = floor(prf * (1/(rpm/60)));
frames = floor(size(p,2)/frameLength);
trim = [1:frames * frameLength];
if(frames==0)
    frames = 1;
    trim = 1:size(p,2);
end

D2 = p(:,trim);
blockTimes2 = blockTimes(:,trim);

D2 = reshape(D2(:),size(D2,1),frameLength,[]);
blockTimes2 = reshape(blockTimes2,frameLength,[]);
%%
pixelSpaceing = 2;
x = 1:1:301;
x = (x-mean(x)) * pixelSpaceing ;

%x=1000;
y =1:1:301;
y = (y-mean(y)) * pixelSpaceing ;
%x = [-50:1:10];
%y =  [-50:1:100];
z = ones(size(y)) ;
clear trialpoint;
count = 1;
%for kk = 0:4:24
for ii = 1:length(x)
    for jj = 1:length(y)

        
        trialpoint(:,count) = [x(ii) y(jj) z(jj)];
        count = count + 1;
        end
  %  end
end
%%
simpletest;
plotsar( P,x,y,80);
