% Test of proper scattering range
% tests that we have the poper range to scatterer

% file parameters
for type = [1]
    addpath('\\BIGSAZENAS\projects\sprdn\jge\ISAR code\')
    addpath('\\BIGSAZENAS\projects\sprdn\jge\ISAR code\test')
    addpath('\\BIGSAZENAS\projects\sprdn\jge\ISAR code\bitmaps')
    
    acount = 1;
    angle = 0;
    rPosition{1}  =20000*[cosd(angle) sind(angle) 1/20];
    rVelocity{1} = 40* [cosd(angle-90) sind(angle-90) 0];
    tPosition{1}=rPosition{1};
    tVelocity{1} = 40* [cosd(angle-90) sind(angle-90) 0];
    angle = 8;
    rPosition{2}  =20000*[cosd(angle) sind(angle) 1/20];
    rVelocity{2} = 40* [cosd(angle-90) sind(angle-90) 0];
    tPosition{2}=rPosition{2};
    tVelocity{2} = 40* [cosd(angle-90) sind(angle-90) 0];
    
    rPosition{3}  =rPosition{1};
    rVelocity{3} = rVelocity{1} ;
    tPosition{3}= tPosition{2};
    tVelocity{3} = rVelocity{2} ;
    for angle = [1:3]
        %%
        
        if(type==0)
            filebase = sprintf('//bigsazenas/data/sprdn/2platformnew/staticsinglesourceresolutiontarget_%d',angle);
        else
            filebase = sprintf('//bigsazenas/data/sprdn/2platformnew/resolutiontargetsource_%d',angle);
        end
        image_filename = 'C:\Users\jeicke\Desktop\sprdn\ISAR code\bitmaps\test3.png';
        %system parameters
        time =8;
        Fs = 900E6; %samples/second
        Fc = 9E9; %hertz
        C = 299792458; %m/s
        
        % RADAR PARAMETERS
        verticalBeamwidth = 25; % in degrees
        horizontalBeamwidth = 25; % in degrees
        prf = 2000; % in hz
        pri = 1/prf;
        pulsSamples = pri * Fs;
        pulseBandwidth = 300E6; % in hz
        dutyFactor = .14;
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
        gainTransmit = 30;
        gainReceive = 30;
        % transmit power in watts
        powerTransmit = 100/dutyFactor;
        
        invert_image = true;
        
        res_outline;
        
        
        
        
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
        
        receiverVelocity = rVelocity{angle};
        receiverPosition =   rPosition{angle};
        receiverPosition = receiverPosition-receiverVelocity * time/2;
        transmitterPosition =tPosition{angle};
        transmitterVelocity = tVelocity{angle};
        transmitterPosition = transmitterPosition-transmitterVelocity * time/2;
        % make transmitter and receiver tracks
        [ track] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );
        
        
        %scatteredpoint = scatteredpoint/9;
        %scatteredpoint(3,:) = scatteredpoint(3,:) * 2;
        
        trackCenter = [0 0 0];
        fprintf(1,'RANGE RESOLUTION %3.2f m\n',C/(2*pulseBandwidth));
        fprintf(1,'CROSS-RANGE RESOLUTION %3.2f m\n',norm(trackCenter-receiverPosition)*lamda/(2*norm(receiverVelocity)*time));
        if(type==0)
            clear scatteredpoint;
            clear reflextivity;
            scatteredpoint(1,1) = 0;
            scatteredpoint(2,1) = 0;
            scatteredpoint(3,1) = 0;
            sigma(1) = 1;
            sigma(2) = 1;
            sigma(3) = 1;
        end
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
        ranges = ([1:1:Parameters.blockSize]-1) * C/Fs * .5+Tz * C * .5-C/Fs*length(pulse_replica)/4+C/Fs/2;
        
        
        
        
        %% Point spread function
        
        xi = [-15:.05:15] + trackCenter(1);
        yi = [-15:.05:15] +  trackCenter(2);
        zi = [0 ] +  trackCenter(3);
        
        [Px] =bpm( filebase,prf,pulse_replica,ranges,track,[0 0 0],xi,yi,zi);
        if(type==0)
            save(sprintf('//bigsazenas/data/sprdn/2platformnew/staticsinglesourceresolutiontarget_%d_sarimage.mat',angle),'Px','xi','yi','zi','prf','pulse_replica','track','Parameters','-v7.3' );
        else
            save(sprintf('//bigsazenas/data/sprdn/2platformnew/resolutiontargetsource_%d_sarimage.mat',angle),'Px','xi','yi','zi','prf','pulse_replica','track','Parameters','-v7.3' );
        end
        %%
        figure
        colormap jet
        imagesc(-xi/1000,yi,20*log10(abs(squeeze(Px))));
        nf = median(20*log10(abs(Px(:))));
        caxis([nf-5 max(20*log10(abs(squeeze(Px(:)))))]);
        ylabel('CROSS-RANGE (M)')
        xlabel('RANGE (KM)');
        hold on;
        colorbar
        plot(squeeze(track(1,1,3:end)),squeeze(track(2,1,3:end)),'ro')
        
        Pxx{acount} = Px;
        acount = acount + 1;
        drawnow
    end
end
return;
%% geometric mean
Pxa = ones(size(Px));
for ii = 1:length(Pxx)
    Pxa = Pxa.*abs(Pxx{ii});
end
Pxa = Pxa.^(1/length(Pxx));
figure
colormap jet
imagesc(-xi/1000,yi,20*log10(abs(squeeze(Pxa))));
nf = median(20*log10(abs(Pxa(:))));
caxis([nf nf+50]);
ylabel('CROSS-RANGE (M)')
xlabel('RANGE (KM)');
hold on;
colorbar
plot(squeeze(track(1,1,3:end)),squeeze(track(2,1,3:end)),'ro')
%% harmonic mean
Pxa = zeros(size(Px));
for ii = 1:length(Pxx)
    Pxa = Pxa+1./abs(Pxx{ii});
end
Pxa = length(Pxx)./Pxa;
figure
colormap jet
imagesc(-xi/1000,yi,20*log10(abs(squeeze(Pxa))));
nf = median(20*log10(abs(Pxa(:))));
caxis([nf nf+50]);
ylabel('CROSS-RANGE (M)')
xlabel('RANGE (KM)');
hold on;
colorbar
plot(squeeze(track(1,1,3:end)),squeeze(track(2,1,3:end)),'ro')
%% mean
Pxa = zeros(size(Px));
for ii = 1:length(Pxx)
    Pxa = Pxa+abs(Pxx{ii});
end
Pxa = Pxa/length(Pxx);
figure
colormap jet
imagesc(-xi/1000,yi,20*log10(abs(squeeze(Pxa))));
nf = median(20*log10(abs(Pxa(:))));
caxis([nf nf+50]);
%caxis([-25 30]);
ylabel('CROSS-RANGE (M)')
xlabel('RANGE (KM)');
hold on;
colorbar
plot(squeeze(track(1,1,3:end)),squeeze(track(2,1,3:end)),'ro')
%% max
Pxa = zeros(size(Px));
for ii = 1:length(Pxx)
    Pxa = Pxa+abs(Pxx{ii});
end

figure
colormap jet
imagesc(-xi/1000,yi,20*log10(abs(squeeze(Pxa/length(Pxx)))));
nf = median(20*log10(abs(Pxa(:)/length(Pxx))));
caxis([nf nf+50]);
%caxis([-25 30]);
ylabel('CROSS-RANGE (M)')
xlabel('RANGE (KM)');
hold on;
colorbar
plot(squeeze(track(1,1,3:end)),squeeze(track(2,1,3:end)),'ro')
%% var
Pxa = zeros(size(Px));

for ii = 1:size(Pxx{1},1)
    for jj = 1:size(Pxx{1},2)
        x = [];
        for kk = 1:length(Pxx)
            x(kk) = Pxx{kk}(ii,jj);
        end
        Pxa(ii,jj) =std(x);
    end
end


figure
colormap jet
imagesc(-xi/1000,yi,20*log10(abs(squeeze(Pxa))));
nf = median(20*log10(abs(Pxa(:))));
caxis([nf nf+50]);
ylabel('CROSS-RANGE (M)')
xlabel('RANGE (KM)');
hold on;
colorbar
plot(squeeze(track(1,1,3:end)),squeeze(track(2,1,3:end)),'ro')
%%

count = 1;
PQ = paddedsize(size(Ppsf{1}));
Ppsf= Pxx;
clear psf B gn;
for ii = 1:length(Pxx)
    psf(:,:,count) = fft2(Ppsf{ii}.^2,PQ(1), PQ(2));
    psf(:,:,count)  = 2*ifft2(psf(:,:,count)/(mean(mean(abs(psf(:,:,count))))));
    B(:,:,count) = abs(Pxx{ii}) ;
    imagesc(20*log10(abs(psf(:,:,count))))
    drawnow
    
    count = count + 1;
    pause(1)
end
g = B;
gn = B;
f = mean(B,3);
fn = mean(B,3);
figure
imagesc(f)
figure
%%
for ii = 1:1000
    z = zeros(size(f));
    for jj = 1:size(psf,3)
        pf = fft2(g(:,:,jj)-gn(:,:,jj),PQ(1), PQ(2));
        pr = ifft(pf.*psf(:,:,jj));
        pr = pr(2:size(Ppsf,1)+1,2:size(Ppsf,2)+1);
        z = z +pr;
    end
    fn = fn+  z;
    imagesc(fn)
    
    for jj = 1:size(psf,3)
        rx = ifft(fft(fn,PQ(1), PQ(2)).*psf(:,:,jj));
        gn(:,:,jj) =  rx(2:size(Ppsf,1)+1,2:size(Ppsf,2)+1);
    end
    
end
%%
figure;
plot(squeeze(track(1,:,1)),squeeze(track(2,:,1)),'ro')
hold on;
plot(squeeze(track(1,:,2)),squeeze(track(2,:,2)),'bo')
plot(squeeze(track(1,:,3)),squeeze(track(2,:,3)),'go')
axis equal
grid on;
xlabel('X (m)')
ylabel('Y (m)')
legend('PLATFORM 1','PLATFORM 2','TARGET')
%%
figure
Px = Pxy{1};

acount = acount + 1;
colormap jet
imagesc(-xi/1000,yi,20*log10(abs(squeeze(Px))));
nf = median(20*log10(abs(Px(:))));
caxis([nf-5 max(20*log10(abs(squeeze(Px(:)))))]);
ylabel('CROSS-RANGE (M)')
xlabel('RANGE (KM)');
hold on;
colorbar
plot(squeeze(track(1,1,3:end)),squeeze(track(2,1,3:end)),'ro')

drawnow