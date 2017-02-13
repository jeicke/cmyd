% Test of proper scattering range
% tests that we have the poper range to scatterer
% RADAR PARAMETERS
addpath('C:\Users\jeicke\Desktop\sprdn\ISAR code\')
addpath('C:\Users\jeicke\Desktop\sprdn\ISAR code\test')
addpath('C:\Users\jeicke\Desktop\sprdn\ISAR code\bitmaps')
USECUDA = false;
Fs = 900E6; %samples/second
Fc = 9E9; %hertz
rpm = 24; %rotations per minute of radar
verticalBeamwidth = 25; % in degrees
horizontalBeamwidth = 25; % in degrees
C = 299792458; %m/s
clear track;
% 200-540
prf = 2000; % in hz
pri = 1/prf;
pulsSamples = pri * Fs;
pulseBandwidth = 600E6; % in hz
dutyFactor = .14;
minimumRange = C*pri*dutyFactor;
pulseLength = pri*dutyFactor ; % in seconds
filebase = 'C:\testdata\test';
image_filename = 'C:\Users\jeicke\Desktop\sprdn\ISAR code\bitmaps\test3.png';
invert_image = true;
%sip_outline;
wife_outline;
gamma = -10;
% noise figure of transmitter in dB
NF = 4;
% losses in dB
losses = 5;
% receiver temperature in Kelvin
T = 290;
% Boltzmann's constant
k  = 1.3806503*10^-23;
gainTransmit = 30;
gainReceive = 30;
% transmit power in watts
powerTransmit = 100/dutyFactor;

lamda = C/Fc;
snr = radarequation(k*T,losses,NF,powerTransmit,gainTransmit,gainReceive,Fs,lamda);

time =12;
% source and receiver locations
receiverVelocity = 50* [0 1 0];
receiverPosition = -receiverVelocity * time/2;
receiverPosition(3) = 0;
transmitterPosition =receiverPosition;
receiverPosition(3) = 0;
transmitterVelocity = receiverVelocity;

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
orientation(3,:) = -90%[0:horizontalBeamwidth:(length(t)-1)*horizontalBeamwidth],360;% - radar_velocity(3)/2  * time;
orientation(4,:) = t;

%geometry of transmit array
[geometry w h]= loadgeometry('rectangle',[],arrayLength,arrayHeight,lamda);

% make transmitter and receiver tracks
[ track] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );

% make scattered point pegged to transmitter (like we're on a boat)

%scatteredpoint = scatteredpoint/9;
%scatteredpoint(3,:) = scatteredpoint(3,:) * 2;

trackCenter = [-15000 0 0];
fprintf(1,'RANGE RESOLUTION %3.2f m\n',C/(2*pulseBandwidth))
fprintf(1,'CROSS-RANGE RESOLUTION %3.2f m\n',norm(trackCenter)*lamda/(2*norm(receiverVelocity)*time))
% clear scatteredpoint;
% clear reflextivity;
% scatteredpoint(1,1) = 0;
% scatteredpoint(2,1) = 0;
% scatteredpoint(3,1) = 0;
% sigma(1) = 1;
% sigma(2) = 1;
% sigma(3) = 1;
for ii = 1:size(scatteredpoint,2);
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
w =hannwindow(length(geometry)/3).';
w = kron(w,[1 1 1]');
w = w/sum(abs(w));
%geometry = [];

%w = [];
%%
Parameters = radsim(filebase,time,track,snr,...
    'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
    'fs',Fs,'fc',Fc,'clutter',reflextivity',...
    'transmitshading',w,'supressdirect',true,'rangecorrect',true );

%%
pulse_replica=m_chirp(Parameters.timeSeriesParameters(1),Parameters.timeSeriesParameters(2),0,Parameters.basebandSamplingRate);

[systemDelay(2), systemDelay(1)] = finddelays( Parameters.track,[0 0 0],Parameters.C,Parameters.timeSeriesParameters(2) );
%[P, blockTimes] =processradarblock( filebase,pulse_replica,systemDelay);
Z = floor((systemDelay(1))*Fs);
Tz = Z/Fs;
ranges = ([1:1:Parameters.blockSize]-1) * C/Fs * .5+Tz * C * .5;




%% Point spread function



windows = {{@rectwin},{@taylorwin,50,-90},{@hann}};
% cross range resolution
for win = 1:length(windows)
    xi = [-20:.075:20] + trackCenter(1);
    yi = [-20:.075:20] +  trackCenter(2);
    zi = [0 ] +  trackCenter(3);
    
    [Px] =bpm( filebase,prf,pulse_replica,ranges,track,[0 0 0],xi,yi,zi,...
        'windowrange',@hann,'windowcrossrange',windows{win},'cuda',true);
    

    figure
    colormap gray
    imagesc(-xi/1000,yi,20*log10(abs(squeeze(Px))));
    nf = median(20*log10(abs(Px(:))));
    caxis([max(caxis)-20 max(caxis)]);
    ylabel('CROSS-RANGE (M)')
    xlabel('RANGE (KM)');
    hold on;
    colorbar
  
    
    
    
end
%%
[Pmt] =bpm( filebase,prf,pulse_replica,ranges,track,[0 0 0],xi,yi,zi,'cuda',false,'multitaper',true,'windowrange',@hann);
%%
 figure
    colormap gray
    imagesc(-xi/1000,yi,20*log10(abs(squeeze(Pmt))));
    nf = median(20*log10(abs(Px(:))));
    caxis([max(caxis)-20 max(caxis)]);
    ylabel('CROSS-RANGE (M)')
    xlabel('RANGE (KM)');
    hold on;
    colorbar

%%
% figure
% colormap jet
% imagesc(-xi/1000,yi,20*log10(abs(squeeze(Px))));
% nf = median(20*log10(abs(Px(:))));
% caxis([max(caxis)-80 max(caxis)]);
% ylabel('CROSS-RANGE (M)')
% xlabel('RANGE (KM)');
% hold on;
% colorbar
    %plot(squeeze(track(1,1,3:end)),squeeze(track(2,1,3:end)),'ro')