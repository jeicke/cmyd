% Test of proper scattering range
% tests that we have the poper range to scatterer

% file parameters

addpath('C:\Users\jeicke\Desktop\sprdn\ISAR code\')
addpath('C:\Users\jeicke\Desktop\sprdn\ISAR code\test')
addpath('C:\Users\jeicke\Desktop\sprdn\ISAR code\bitmaps')
acount = 1;
for angle = [0:20:360]
    filebase = 'C:\testdata\test';
    image_filename = 'C:\Users\jeicke\Desktop\sprdn\ISAR code\bitmaps\test3.png';
    %system parameters
    time =2;
    Fs = 200E3; %samples/second
    Fc = 1.3E9; %hertz
    C = 299792458; %m/s
    
    % RADAR PARAMETERS
    verticalBeamwidth = 25; % in degrees
    horizontalBeamwidth = 25; % in degrees
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
    sip_outline;
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
    
    receiverVelocity = 60* [cosd(angle-90) sind(angle-90) 0];
    receiverPosition  =[cosd(angle) sind(angle) 0];
    receiverPosition = receiverPosition-receiverVelocity * time/2;
    receiverPosition(3) = 10E3;
    transmitterPosition =[-32.5; 0; 20]*1000;

    transmitterVelocity = [0 0 0];
    
    % make transmitter and receiver tracks
    [ track] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );
    
    
    %scatteredpoint = scatteredpoint/9;
    %scatteredpoint(3,:) = scatteredpoint(3,:) * 2;
    
    trackCenter = [0 0 0];
    fprintf(1,'RANGE RESOLUTION %3.2f m\n',C/(2*pulseBandwidth))
    fprintf(1,'CROSS-RANGE RESOLUTION %3.2f m\n',norm(trackCenter-receiverPosition)*lamda/(2*norm(receiverVelocity)*time))
    %     clear scatteredpoint;
    %     clear reflextivity;
    %     scatteredpoint(1,1) = 0;
    %     scatteredpoint(2,1) = 0;
    %     scatteredpoint(3,1) = 0;
    %     sigma(1) = 1;
    %     sigma(2) = 1;
    %     sigma(3) = 1;
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
    %     w =hannwindow(length(geometry)/3).';
    %     w = kron(w,[1 1 1]');
    %     w = w/sum(abs(w));
    %geometry = [];
    
    %w = [];
    %%
    Parameters = radsim(filebase,time,track,snr,...
        'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
        'fs',Fs,'fc',Fc,'clutter',reflextivity','supressdirect',true,'rangecorrect',false );
    
    %%
    pulse_replica=m_chirp(Parameters.timeSeriesParameters(1),Parameters.timeSeriesParameters(2),0,Parameters.basebandSamplingRate);
    
    [systemDelay(2), systemDelay(1)] = finddelays( Parameters.track,[0 0 0],Parameters.C,Parameters.timeSeriesParameters(2) );
    %[P, blockTimes] =processradarblock( filebase,pulse_replica,systemDelay);
    Z = floor((systemDelay(1))*Fs);
    Tz = Z/Fs;
    ranges = ([1:1:Parameters.blockSize]-1) * C/Fs * .5+Tz * C * .5;
    
    
    
    
    %% Point spread function
    
    xi = [-55:.5:55] + trackCenter(1);
    yi = [-50:.5:50] +  trackCenter(2);
    zi = [0 ] +  trackCenter(3);
    
    [Px] =bpm( filebase,prf,pulse_replica,ranges,track,[0 0 0],xi,yi,zi,...
        'windowrange',@hann,'windowcrossrange',@hann);
    
    
    figure
    colormap jet
    imagesc(-xi/1000,yi,20*log10(abs(squeeze(Px))));
    nf = median(20*log10(abs(Px(:))));
    caxis([nf-5 nf+50]);
    ylabel('CROSS-RANGE (M)')
    xlabel('RANGE (KM)');
    hold on;
    colorbar
    plot(squeeze(track(1,1,3:end)),squeeze(track(2,1,3:end)),'ro')
    
    Pxx{acount} = Px;
    acount = acount + 1;
    
end
%%
Pxa = zeros(size(Px));
for ii = 1:length(Pxx)
    Pxa = Pxa+abs(Pxx{ii});
end
figure
colormap jet
imagesc(-xi/1000,yi,20*log10(abs(squeeze(Pxa/length(Pxx)))));
nf = median(20*log10(abs(Pxa(:)/length(Pxx))));
%caxis([-25 30]);
ylabel('CROSS-RANGE (M)')
xlabel('RANGE (KM)');
hold on;
colorbar
plot(squeeze(track(1,1,3:end)),squeeze(track(2,1,3:end)),'ro')

