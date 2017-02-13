% Test of proper scattering range
% tests that we have the poper range to scatterer
% RADAR PARAMETERS

Fs = 60E6; %samples/second
Fc = 9410E6; %hertz
rpm = 24; %rotations per minute of radar
verticalBeamwidth = 25; % in degrees
horizontalBeamwidth = 5; % in degrees
prf = 2404; % in hz
pulseBandwidth = 30E6; % in hz
pulseLength = 4800E-9; % in seconds
C = 299792458; %m/s
filebase = 'C:\Users\rnl\Desktop\testdata\test';
image_filename = 'C:\Users\rnl\Desktop\ISAR code\bitmaps\test5.png';
invert_image = true;
shipmask;
time =1.2/(rpm/60);
% source and receiver locations
receiverVelocity = 200* [0 1 0];   
receiverPosition = -receiverVelocity * time/2; 
receiverPosition(3) = 1000;
transmitterPosition = [10000 0 45];%receiverPosition ;%[10000 0 0];

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
w =hanning(length(geometry)/3);
w = kron(w,[1 1 1]');
w = w/sum(abs(w));
%geometry = [];
% w = [];
%%
radsim(filebase,time,track,0,...
       'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
       'fs',Fs,'fc',Fc,'transmitgeometry',geometry,...
       'transmitorientation',orientation,'clutter',reflextivity',...
       'transmitshading',w); 
%%
pulse_replica=m_chirp(pulseBandwidth,pulseLength,0,Fs)'.';
 
pulse_replica=pulse_replica(end:-1:1);
rngspa= 299792458/Fs;
ranges = [0:1:plength] * rngspa- (length(pulse_replica)-2) * rngspa;

rm = norm(receiverPosition -transmitterPosition);
rm = floor(rm/rngspa)+length(pulse_replica)-1;

step = floor(1500 / rngspa);

 
[maxDelay minDelay ] = finddelays( track,geometry,C );

[h minRangeIndex ] = min(abs(ranges-minDelay *C));
% minRangeIndex  = 3000;
[h maxRangeIndex ] = min(abs(ranges-maxDelay *C));
 
rangeGates = (max([1  minRangeIndex-2*length(pulse_replica)]):min([maxRangeIndex+2*length(pulse_replica) length(ranges)]));
%rangeGates = [max([1 rm-step]):1:min([rm+step length(ranges)])];
%rangeGates = 1:2000;
[D blockTimes] = processradar( filebase,prf,1,pulse_replica,rangeGates );
rangeShift = ranges(1)-ranges(rangeGates(1));
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
x = 1:1:201;
x = (x-mean(x)) * pixelSpaceing ;
x =  [-100:1:100];
%x=1000;
y =1:1:201;
y = (y-mean(y)) * pixelSpaceing ;
x = [-100:3:450];
y =  [-100:3:100];
z = ones(size(y)) *10;
clear trialpoint;
count = 1;
for ii = 1:length(x)
    for jj = 1:length(y)
        
        
        trialpoint(:,count) = [x(ii) y(jj) z(jj)];
        count = count + 1;
    end
end

