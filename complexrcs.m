% Test of proper scattering range
% tests that we have the poper range to scatterer
% RADAR PARAMETERS

Fs = 20E6; %samples/second
Fc = 9410E6; %hertz
rpm = 24; %rotations per minute of radar
verticalBeamwidth = 25; % in degrees
horizontalBeamwidth = 1.8; % in degrees
prf = 2304; % in hz
pulseBandwidth = 18E6; % in hz

C = 299792458; %m/s
filebase = 'C:\Users\rnl\Desktop\testdata\test';
invert_image = true;

% scatteredpoint = [];
% scatteredpoint(:,1) = [20000 0 0];
% scatteredpoint(:,2) = [40000 0 0];
% sigma= [];
% sigma(1) = 1;
% sigma(2) = 1;
% constant gamma cluter
gamma = -10;
% noise figure of transmitter in dB
NF = 3;
% losses in dB
losses = 5;
% receiver temperature in Kelvin
T = 290;
% Boltzmann's constant
k  = 1.3806503*10^-23;
gainTransmit = 28;
gainReceive = 0;
% transmit power in watts
powerTransmit = 200;
%dutyCycle in percent
dutyCycle = 10;
pulseLength =1/prf  * dutyCycle * .01 ; % in seconds
lamda = C/Fc;
snr = radarequation(k*T,losses,NF,powerTransmit,gainTransmit,gainReceive,Fs,lamda);
shipmask;
sigma = sigma;

time =1.5/(rpm/60);
% source and receiver locations
receiverVelocity =200* [0 1 0];   
receiverPosition = -receiverVelocity * time/2; 
receiverPosition(3) = 0;
transmitterPosition = [1000 0 45];%receiverPosition ;%[10000 0 0];

transmitterVelocity=[0 0 0];%receiverVelocity;%[0 0 0];

% derived parameters

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
[geometryVertical w h]= loadgeometry('line',[],arrayHeight*2,1,lamda);
geometryVertical(:,3) = geometryVertical(:,1);
geometryVertical(:,1) = 0;
[geometryHorizontal w h]= loadgeometry('line',[],arrayLength,1,lamda);
geom{2} = geometryVertical;
geom{1} = geometryHorizontal;

% make transmitter and receiver tracks
[ track] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );
[ track2] =maketrack(time,time/10,0*transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );

% make scattered point pegged to transmitter (like we're on a boat)



clear reflextivity;
for ii = 1:size(scatteredpoint,2);
    %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    
    reflextivity(ii) = sigma(ii);
    track(:,:,ii+2) = track(:,:,1);
   % if(ii>200)
        %track(:,:,ii+2) = track2(:,:,1);
   % end
    track(1,:,ii+2) =track(1,:,ii+2) +scatteredpoint(1,ii);
    track(2,:,ii+2) =track(2,:,ii+2) +scatteredpoint(2,ii);
    track(3,:,ii+2) =track(3,:,ii+2) +scatteredpoint(3,ii);
    
end
snr = 10*log10(snr) *ones(1,length(sigma)+1);
w =hanning(length(geometryHorizontal));
w = w/sum(abs(w));

weights{1} = w;
w =chebwin(size(geometryVertical,1));
w = w/sum(abs(w));
weights{2} = w;

[geometry w h]= loadgeometry('rectangle',[],arrayLength,arrayHeight,lamda);
w =hanning(length(geometry)/3);
w = kron(w,[1 1 1]');
w = w/sum(abs(w));

w = kron(w,chebwin(3));

%geom = [];
 %w = [];
%%
radsim(filebase,time,track,snr,...
       'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
       'fs',Fs,'fc',Fc,'transmitgeometry',geom,...
       'transmitorientation',orientation,'clutter',reflextivity',...
       'transmitshading',weights,'rangecorrect',false); 
%%
[maxDelay minDelay ] = finddelays( track,geometryVertical,C );
pulse_replica=m_chirp(pulseBandwidth,pulseLength,0,Fs)'.';
 
pulse_replica=pulse_replica(end:-1:1);
rngspa= 299792458/Fs;
ranges = [0:1:(plength-1)] * rngspa- (length(pulse_replica)-2) * rngspa;
drange = ranges(end)-ranges(1);
while ranges(end)< (minDelay*C)
    ranges(end)
    ranges = ranges+drange;
end  
rm = norm(receiverPosition -transmitterPosition);
rm = floor(rm/rngspa)+length(pulse_replica)-1;

step = 1000;%floor(1500 / rngspa);

 


[h minRangeIndex ] = min(abs(ranges-minDelay *C));
% minRangeIndex  = 3000;
[h maxRangeIndex ] = min(abs(ranges-maxDelay *C));
 
rangeGates = (max([1  minRangeIndex-length(pulse_replica)]):min([maxRangeIndex+length(pulse_replica) length(ranges)]));
rangeGates = [max([1 rm-step]):1:min([rm+step length(ranges)])];
%rangeGates = 1:2000;
[D blockTimes] = processradar( filebase,prf,1,pulse_replica,rangeGates,[],1  );
f = (prf/2)* (-size(D,2)/2:size(D,2)/2-1)/(size(D,2)/2);

%[D ] = pulsecompression(D,pulse_replica,'cheb',70);
dt = 1/Fs;
% p1 = 10*log10(dt*trapz(abs(D(:,:,1)).^2))
% 
% 
% p2 = 10*log10(dt*trapz(abs(D2(:,:,1)).^2)/2)
%integrate D and D2 over pulse

rangeShift = ranges(1)-ranges(rangeGates(1));
ranges = ranges(rangeGates );


% plot beampattern
%%
plottrack(track,[]);
figure;
p = 20 * log10(abs(squeeze(D(860:1000,:,:))));
imagesc(aztransmit(50,:) * 180/pi,ranges(860:1000),p)
ylabel('RANGE (M)')
xlabel('BEARING (DEG)')
caxis([max(s(:))-60 max(s(:))])
colorbar
figure;
[r_transmitter_scatter aztransmit detransmit]= computeRange(track,  blockTimes(1,:),2);
 bp = sum(abs(D(950:960,:)),1);


[beamPower] = beampattern(Fc,aztransmit(50,:) * 180/pi,...
                zeros(size(aztransmit(50,:))),ones(size(scatteredpoint,2),1)',[], [],scatteredpoint',[],C);
  bp1 = 20*log10(bp);
  bp2 = 10*log10(beamPower);
  plot(aztransmit(50,:) * 180/pi,bp1-max(bp1),'r',aztransmit(50,:) * 180/pi,bp2-max(bp2),'linewidth',1.5)
  xlabel('ANGLE (DEG)')
  ylabel('POWER (dB)')
  grid
  axis([min(aztransmit(50,:) * 180/pi) max(aztransmit(50,:) * 180/pi) -50 5])
  legend('SCATTER MODEL','THEORY')
%%
threeD = false;
if(~threeD)
    p = squeeze(sum(D,2));
   % p = squeeze(D(:,17,:));
    dt = 1/prf;
    t = ([0:1:(size(D,3)-1)])*(dt * size(D,2))+dt/2;
    azx = interp1(orientation(4,:),orientation(3,:),t);
    
    
    
    frameLength = floor(prf * (1/(rpm/60))/size(D,2));
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
else
    [D1] = dopplerfilter(D,32 ,'cheb',70,false);
    dt = 1/prf;
    t = ([0:1:(size(D1,3)-1)])*(dt * size(D1,2))+dt/2;
    azx = interp1(orientation(4,:),orientation(3,:),t);
    
    
    
    frameLength = floor(prf * (1/(rpm/60))/size(D1,2));
    frames = floor(size(D1,3)/frameLength);
    trim = [1:frames * frameLength];
    if(frames==0)
        frames = 1;
        trim = 1:size(D1,3);
    end
    
    D2 = D1(:,:,trim);
    %D2 = reshape(D2(:),size(D2,1),frameLength,[]);
    %blockTimes2 = reshape(blockTimes2,frameLength,[]);
end
%%
pixelSpaceing = 2;
x = 1:1:201;
x = (x-mean(x)) * pixelSpaceing ;
x =  [-100:1:100];
%x=1000;
y =1:1:201;
y = (y-mean(y)) * pixelSpaceing ;
x = [-100:1:500];
y =  [-100:2:100];
z = ones(size(y)) *10;
clear trialpoint;
count = 1;
for ii = 1:length(x)
    for jj = 1:length(y)
        
        
        trialpoint(:,count) = [x(ii) y(jj) z(jj)];
        count = count + 1;
    end
end



