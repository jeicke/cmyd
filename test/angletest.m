% Test of proper scattering range
% tests that we have the poper range to scatterer

% file parameters

acount = 1;

filebase = '\\hermes\data\cmdn\testdata';
%image_filename = 'C:\Users\jeicke\Desktop\sprdn\ISAR code\bitmaps\test3.png';
%system parameters
time =2;
Fs = 600E3; %samples/second
Fc = 1.3E9; %hertz
C = 299792458; %m/s

% RADAR PARAMETERS

prf = 1000; % in hz
pri = 1/prf;
pulsSamples = pri * Fs;
pulseBandwidth = 100E3; % in hz
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



%boat_outline;

lamda = c/operating_frequency ;
% width of array in meters
width = 6*lamda;
% height of array in meters
height = lamda;
% constant gamma cluter
gamma = -10;
% 3:1 overlap (75%)
overlap = 0;

[fullarrayGeometry]= load_geometry('line',[],width,height,299792458/operating_frequency );
%fullarrayGeometry = [0 0 0];
lamda = C/Fc;
snr = radarequation(k*T,losses,NF,powerTransmit,gainTransmit,gainReceive,Fs,lamda);
% derived parameters
lamda = C/Fc;
arrayHeight = (1.27 * lamda)/(verticalBeamwidth *pi/180);
arrayLength = (1.27 * lamda)/(horizontalBeamwidth *pi/180);
% dwell_time = 1/(rpm/60) * horizontalBeamwidth/360;
Npulses = floor(time * prf);%number of pulse in look direction
plength = round(Fs/ prf);
clear track;

% receiver and transmitter locations
angle = 90;
receiverVelocity = 0* [cosd(angle-90) sind(angle-90) 0];
receiverPosition  =[cosd(angle) sind(angle) 0];
receiverPosition = receiverPosition-receiverVelocity * time/2;
receiverPosition(3) = 10E3;
transmitterPosition =[-20 0 10]*1000;

transmitterVelocity =[0 0 0];

% make transmitter and receiver tracks
[ track] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );


%scatteredpoint = scatteredpoint/9;
%scatteredpoint(3,:) = scatteredpoint(3,:) * 2;
ANG = -25;
trackCenter = [3000*sind(ANG) 3000*cosd(ANG) 10000];
fprintf(1,'RANGE RESOLUTION %3.2f m\n',C/(2*pulseBandwidth))
fprintf(1,'CROSS-RANGE RESOLUTION %3.2f m\n',norm(trackCenter-receiverPosition)*lamda/(2*norm(receiverVelocity)*time))
clear scatteredpoint;
clear reflextivity;
scatteredpoint(1,1) = 0;
scatteredpoint(2,1) = 0;
scatteredpoint(3,1) = 0;
sigma(1) = 1;
sigma(2) = 1;
sigma(3) = 1;
for ii = 1:size(scatteredpoint,2)
    %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    
    reflextivity(ii) = sigma(ii);
    track(:,:,ii+2) = track(:,:,1);
    track(1,:,ii+2) =trackCenter(1);
    track(2,:,ii+2) =trackCenter(2);
    track(3,:,ii+2) =trackCenter(3);
    %track(1:3,:,ii+2) = 0;
    track(1,:,ii+2) =track(1,:,ii+2)+scatteredpoint(1,ii);
    track(2,:,ii+2) =track(2,:,ii+2)+scatteredpoint(2,ii);
    track(3,:,ii+2) =track(3,:,ii+2)+scatteredpoint(3,ii);
    
end

%reflextivity = reflextivity(1:size(track,3)-1);
snr = 10*log10(snr) *ones(1,length(sigma)+1);


Parameters = radsim(filebase,time,track,snr,...
    'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
    'fs',Fs,'fc',Fc,'clutter',reflextivity','supressdirect',true,'rangecorrect',false,'geometry',fullarrayGeometry );

  pulse_replica=m_chirp(Parameters.timeSeriesParameters(1),Parameters.timeSeriesParameters(2),0,Parameters.basebandSamplingRate);
[systemDelay(2), systemDelay(1)] = finddelays( Parameters.track,[0 0 0],Parameters.C,Parameters.timeSeriesParameters(2) );
%[P, blockTimes] =processradarblock( filebase,pulse_replica,systemDelay);
Z = floor((systemDelay(1))*Fs);
Tz = Z/Fs;
ranges = ([1:1:Parameters.blockSize]-1) * C/Fs * .5+Tz * C * .5;


 [P] =stap( filebase,prf,32,pulse_replica,ranges,track,fullarrayGeometry ,xi,yi,zi);


for ii = 1
    ii
    hold off;
    colormap jet;
    imagesc(fftshift(20*log10(abs(P(:,:,end,ii).')),2));
    caxis([20 60]);
    colorbar;
    drawnow; 
end


az = -90:90;
el = zeros(size(az));
x = squeeze(P(1,540,:,1));
[powerD] = beampattern(Fc,az,el,x,[],[],fullarrayGeometry,[0 0],C);
[powerT] = beampattern(Fc,az,el,[],[],[],fullarrayGeometry,[ANG 0],C);
Pz = 20*log10(abs(powerD));
Pz = Pz-max(Pz);
close all;
plot(az,Pz);
hold on;
plot(az,20*log10(abs(powerT)),'b--');
return
%%
close all
 plottrack(track,[],false,true)
 axis equal
