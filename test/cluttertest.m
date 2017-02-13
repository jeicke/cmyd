% initial test of clutter model
% write simulation files to this directort
filebase = '\\hermes\data\cmdn\testdata';
%time in simulation
time =1;
%load radar parameters
radarparameters;
% make track of transmitter and receiver
createtrack;
%make clutter model
cluttermodel;
%make simulated time series
%%
Parameters = radsim(filebase,time,track,1,...
    'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
    'fs',Fs,'fc',Fc,'clutter',reflextivity','supressdirect',true,'rangecorrect',true,'geometry',fullarrayGeometry );
%make pulse reference
pulse_replica=m_chirp(Parameters.timeSeriesParameters(1),Parameters.timeSeriesParameters(2),0,Parameters.basebandSamplingRate);
[systemDelay(2), systemDelay(1)] = finddelays( Parameters.track,[0 0 0],Parameters.C,Parameters.timeSeriesParameters(2) );
%[P, blockTimes] =processradarblock( filebase,pulse_replica,systemDelay);
Z = floor((systemDelay(1))*Fs);
Tz = Z/Fs;
ranges = ([1:1:Parameters.blockSize]-1) * C/Fs * .5+Tz * C * .5;
%% do stap

[P] =stap( filebase,prf,32,pulse_replica,ranges,track,fullarrayGeometry ,xi,yi,zi,'windowrange',{@chebwin,40},'windowdoppler',{@chebwin,90});
D = P(:,250:330,:,ii);
D = permute(D,[2 1 3]);
angle = 0;
Vs = exp(-1i*2*pi*1/lamda * fullarraGeometry*[cosd(angle);sind(angle);0]);
Vd = dftmtx(32);

[Dbf, w] = beamform_stap(D,Vs,Vd,20,1,.02,5);
 %[P] =stap( filebase,prf,256,pulse_replica,ranges,track,fullarrayGeometry ,xi,yi,zi,'windowdoppler',{@chebwin,90});
%% plot results
for angle = [-90:5:90]
    
    sv = exp(-1i*2*pi*1/lamda * fullarraGeometry*[cosd(angle);sind(angle);0]);
    for ii = 1%:size(P,4)
        ii
        hold off;
        colormap jet;
        Px = P(:,250:330,:,ii);
        Px = reshape(Px,size(Px,1)*size(Px,2),size(Px,3))*sv;
        Px = reshape(Px,size(P,1),[]);
        z = fftshift(20*log10(abs(sum(Px ,3).')),2);
        figure(1)
        imagesc(z);
        colorbar;
       
        figure(2)
        plot(z(20,:))
        % caxis([20 90]);
         drawnow;
    end
end
%%
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
%% plot results
for angle = AZ
    sv = exp(1i*2*pi*1/lamda * fullarraGeometry*[sind(angle);cosd(angle);0]);
    angle
    for ii = 1
        ii
        hold off;
        colormap jet;
        Px = P(:,520:600,1,ii);
%         Px = reshape(Px,size(Px,1)*size(Px,2),size(Px,3));
%         Pn = reshape(Px,size(P,1),[]);
        z = fftshift(20*log10(abs(Px)),2);
        imagesc(20*log10(abs(Px)));
        % caxis([20 90]);
        colorbar;
        drawnow;
    end
end