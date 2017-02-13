% initial test of clutter model
% write simulation files to this directort
filebase = 'C:\Users\jeicke\Documents\data\';
%time in simulation
time =.25;
%load radar parameters
radarparameters;
% make track of transmitter and receiver
createtrack;
%make clutter model
cluttermodel;
%make target
%maketarget;
%create theoretical sinr
%%
close all
 plottrack(track,[],false,true)
 hold on;
 plot(fullarrayGeometry(:,1)*100,fullarrayGeometry(:,2)*100,'k')
 axis equal
%%
Parameters = radsim(filebase,time,track,1,...
    'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
    'fs',Fs,'fc',Fc,'clutter',reflextivity','supressdirect',true,'rangecorrect',true,'geometry',fullarrayGeometry );
%make pulse reference
%%
close all;
[sinr,dopplerbins1,Rx]= theorysinr(trackt,fullarrayGeometry ,prf,32,lamda,reflextivityt,transmitterVelocity,receiverVelocity);
vc = 1;
%vc =lamda/2;
hold on;
for ii = 1:4%length(sinr)
    plot(vc*dopplerbins1,10*log10(abs(sinr{ii})),'linewidth',2)
end
grid on
ylabel('SINR (dB)');
xlabel('SPEED (m/s)')
%%
az = 0:180;
el = zeros(size(az));
[powerD] = beampattern(Fc,az,el,[],[],[],fullarrayGeometry,[90 0],C);
figure;
plot(az,20*log10(abs(powerD)));
%legend('FULLY ADAPTIVE','CONVENTIONAL','POST-DOPPLER');
%make simulated time series

%%
pulse_replica=m_chirp(Parameters.timeSeriesParameters(1),Parameters.timeSeriesParameters(2),0,Parameters.basebandSamplingRate);
[systemDelay(2), systemDelay(1)] = finddelays( Parameters.track,[0 0 0],Parameters.C,Parameters.timeSeriesParameters(2) );
%[P, blockTimes] =processradarblock( filebase,pulse_replica,systemDelay);
Z = floor((systemDelay(1))*Fs);
Tz = Z/Fs;
ranges2 = ([1:1:Parameters.blockSize]-1) * C/Fs +Tz * C-C*length(pulse_replica)/Fs;
%% do stap
vc =1;% lamda;
Vs = dftmtx(size(fullarrayGeometry,1));
Vd = dftmtx(64);
[P,sinrEx,dopplerbins] =stap( filebase,prf,128,pulse_replica,ranges2,track,fullarrayGeometry ,Vs,Vd,30,1,receiverVelocity,Fc, 'windowrange',{@chebwin,40});
[P,sinrEx2,dopplerbins] =stap( filebase,prf,128,pulse_replica,ranges2,track,fullarrayGeometry ,Vs,Vd,30,1,receiverVelocity,Fc,'dopplerwarping',ranges);
figure(1)

hold off;
plot(vc*dopplerbins,10*log10(abs(sinrEx{1})),'b','linewidth',3);
hold on;
plot(vc*dopplerbins,10*log10(abs(sinrEx{2})),'k','linewidth',3);
plot(vc*dopplerbins,10*log10(abs(sinrEx{3})),'r','linewidth',3);
plot(vc*dopplerbins,10*log10(abs(sinrEx{4})),'c','linewidth',3);
grid on
ylabel('SINR (dB)');
xlabel('SPEED (m/s)')
figure(2)
hold off;
plot(vc*dopplerbins,10*log10(abs(sinrEx2{1})),'b','linewidth',3)
hold on;
plot(vc*dopplerbins,10*log10(abs(sinrEx2{2})),'k','linewidth',3)
plot(vc*dopplerbins,10*log10(abs(sinrEx2{3})),'r','linewidth',3)
plot(vc*dopplerbins,10*log10(abs(sinrEx2{4})),'c','linewidth',3)
grid on
ylabel('SINR (dB)');
xlabel('SPEED (m/s)')
 %[P] =stap( filebase,prf,256,pulse_replica,ranges,track,fullarrayGeometry ,xi,yi,zi,'windowdoppler',{@chebwin,90});

%%
pulseInCpi = 128;
vc =1;% lamda;
doppler_bins =vc*linspace(-prf/2,prf/2,pulseInCpi);
[Pn] =rps( filebase,prf,pulseInCpi ,pulse_replica,ranges,track,fullarrayGeometry, 'windowrange',{@chebwin,40},'windowdoppler',{@chebwin,80});
Pn = Pn{1};
%for angle = [-90:10:90]
    angle = 90;
    sv = exp(-1i*2*pi*1/lamda * fullarrayGeometry*[sind(angle);cosd(angle);0]);
    for ii = 1:size(Pn,4)
        ii
        hold off;
        colormap jet;
        Px = Pn(:,:,:,ii);
        Px = reshape(Px,size(Px,1)*size(Px,2),size(Px,3))*sv;
        Px = reshape(Px,size(Pn,1),[]);
       zc = fftshift(20*log10(abs(sum(Px ,3).')),2);
        %zc = 20*log10(abs(sum(Px ,3).'));
        figure(1)
        imagesc(doppler_bins,[1 12.5],zc);
         caxis([-160 -40]);
        colorbar;
       
%         figure(2)
%         hold on;
%         plot(zc(290,:))
%         
%          drawnow;s
%            figure(3)
%            hold on;
   % plot(zc(:,4))
    % caxis([20 90]);
    drawnow
    end
    %%
    [r_scatter_ant, azreceive, dereceive]= computeRange(track, 0*pulseInCpi *1/prf,2);
    figure
    az = -90:90;
    el = zeros(size(az));
    for ii = 280:300
        x = squeeze(Pn(4,ii,:,1));
        w = x;
        ANG = azreceive(end)*180/pi;
        de = dereceive(end) * 180/pi;
        [powerD] = beampattern(Fc,az,el,x,[],[],fullarrayGeometry,[0 0],C);
        [powerW] = beampattern(Fc,az,el,w,[],[],fullarrayGeometry,[0 0],C);
        
        [powerT] = beampattern(Fc,az,el,[],[],[],fullarrayGeometry,[ANG de],C);
        Pz = 20*log10(abs(powerD));
        Pz = Pz-max(Pz);
        hold off
        plot(az,Pz);
        hold on;
        plot(az,20*log10(abs(powerT)),'b--');
        drawnow
        
    end
%end
%% now staggered pri postDoppler
snapshots = 20;
guardband = 1;
sigma = 100;
[Pn] =rps( filebase,prf,pulseInCpi ,pulse_replica,ranges,track,fullarrayGeometry, 'windowrange',{@chebwin,40},'windowdoppler',{@chebwin,80});
[Pf] =rps( filebase,prf,pulseInCpi ,pulse_replica,ranges,track,fullarrayGeometry, 'windowrange',...
    {@chebwin,40},'windowdoppler',{@chebwin,80},'nofft',true);
Pn = Pn{1};
Pf = Pf{1};
clear P;
wd = chebwin(pulseInCpi,80);
wd = wd./norm(wd);
for angle = [90]
    %spatial steering vector

    clear sinr;
    sv = exp(1i*2*pi*1/lamda * fullarrayGeometry*[sind(angle);cosd(angle);0]);
    Vd = dftmtx(pulseInCpi );
    BOM = kron(sv ,Vd);
    for ii = 1%:size(Pn,4)
        
        hold off;
        colormap jet;
        Px = Pn(:,:,:,ii);
        Pxf = Pf(:,:,:,ii);
        count = 1;
        
        for rg = 200:400%1:size(Px,2);
            for dpbin = 1:size(Px,1)
                
                % estimate spatial-doppler correlation matrix for each range using N snapshots
                % around given range value with a given guard band
                startLeft = max([rg - (snapshots+guardband) 1 ]);
                endLeft = max([rg  - guardband-1 1]);
                startRight = min([rg  + guardband+1 size(Px,2)-1]);
                endRight = min([rg +(snapshots+guardband) size(Px,2)-1]);
                data =Pxf(rg,:);
                s = squeeze((Px(dpbin,[startLeft:endLeft startRight:endRight],:)))'.';
                R = (s' * s)/size(s,1);
                R = R + mean(diag(R))*sigma*eye(size(R));
                w = inv(R)*sv;
                w = w./dot(w,R*w);
                wa = kron(w'.',Vd(:,dpbin).*wd);
                v = squeeze((Px(dpbin,rg,:)));
                P(count,dpbin) = abs(w'*v);
                P2(count,dpbin) = abs(wa'*data');
                sinr(count,dpbin) = abs(wa'*BOM(:,dpbin)).^2./abs(wa'*data').^2;
            end
            count = count + 1;
        end
   
    figure(4)
    colormap jet
    imagesc(doppler_bins,[1 12.5],12.5+fftshift(20*log10(abs(P)),2))
    colorbar
    %caxis([-120 -40]);
    drawnow;
    figure(5)
    colormap jet
    imagesc(doppler_bins,[1 12.5],12.5+fftshift(20*log10(abs(P2)),2))
    colorbar
    %caxis([-120 -40]);
    drawnow;
end
end
%%
    [r_scatter_ant, azreceive, dereceive]= computeRange(track, 0*pulseInCpi *1/prf,2);
    figure
    az = -90:.1:90;
    el = zeros(size(az));
    for ii = 280
        x = squeeze(Pn(2,ii,:,1));
        
        ANG = azreceive(end)*180/pi;
        [powerD] = beampattern(Fc,az,el,x,[],[],fullarrayGeometry,[0 0],C);
        [powerW] = beampattern(Fc,az,el,w,[],[],fullarrayGeometry,[0 0],C);
        
        [powerT] = beampattern(Fc,az,el,[],[],[],fullarrayGeometry,[ANG 0],C);
        Pz = 20*log10(abs(powerD));
        Pz = Pz-max(Pz);
        hold off
        plot(az,Pz);
        Pw = 20*log10(abs(powerW));
        Pw = Pw-max(Pw);
        hold on;
        plot(az,Pw,'r');
        hold on;
        %plot(az,20*log10(abs(powerT)),'b--');
        drawnow
        
    end
%%

snapshots = 10;
guardband = 5;
sigma = 10000;
staggers = 1;
[Pn] =rps( filebase,prf,pulseInCpi ,pulse_replica,ranges,track,fullarrayGeometry, 'windowrange',{@chebwin,40},'windowdoppler',{@chebwin,80},'postdoppler',staggers);

% make temporal steering vector matrix
N = size(Pn{1},1);
dsv = dftmtx(N);
clear P Px;
for angle = [0]
    %spatial steering vector
    sv = exp(1i*2*pi*1/lamda * fullarrayGeometry*[cosd(angle);sind(angle);0]);
    Vd = dftmtx(pulseInCpi );
    BOM = kron(sv ,Vd);
    
    In = eye(size(fullarrayGeometry,1));
    clear P;
    for ii = 1:size(Pn{1},4)
        [r_scatter_ant, azreceive, dereceive]= computeRange(track, (ii-1)*pulseInCpi *1/prf,2);
        azreceive(end)*180/pi
        sv = exp(1i*2*pi*1/lamda * fullarrayGeometry*[cosd(azreceive(end));sind(azreceive(end));0]);
        Vd = dftmtx(pulseInCpi );
    BOM = kron(sv ,Vd);
    
    In = eye(size(fullarrayGeometry,1));
        ii
        hold off;
        colormap jet;
        for jj = 1:numel(Pn)
            Px{jj} = Pn{jj}(:,:,:,ii);
        end
        for rg = 200:400;%size(Px,2);
            
            for dpbin = 1:(size(Px{1},1)-staggers +1)
                
                BOM = kron(sv ,Vd(:,dpbin ));
                % estimate spatial-doppler correlation matrix for each range using N snapshots
                % around given range value with a given guard band
                startLeft = max([rg - (snapshots+guardband) 1 ]);
                endLeft = max([rg  - guardband-1 1]);
                startRight = min([rg  + guardband+1 size(Px,2)-1]);
                endRight = min([rg +(snapshots+guardband) size(Px,2)-1]);
                s = [];
                
                Fm = zeros(length(dsv),numel(Pn));
                for pdb =  1:numel(Pn)
                    s = [s squeeze((Px{pdb}(dpbin,[startLeft:endLeft startRight:endRight],:)))];
                    endIndex = pulseInCpi -(numel(Pn)-pdb);
                    Fm(pdb:endIndex,pdb) = dsv(:,dpbin);
                end
                T = kron(Fm,In);
                
                R = (s' * s)/size(s,1);
                R = R + mean(diag(R))*sigma*eye(size(R));
                
                w = inv(R)*(T'*BOM);
               
                w = w./dot(w,R*w);
                 
                 v = [];
                 for pdb =  1:numel(Pn)
                     v = [v squeeze((Px{pdb}(dpbin,rg,:))).'];
                 end
                 v = v.';
                 
                P(rg-199,dpbin) = abs(w'*v);
            end
        end
   
    figure(4)
    colormap jet
    imagesc(12.5+fftshift(20*log10(abs(P)),2))
    colorbar
    caxis([-120 -40]);
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
pulseInCpi = 32;
[Pn] =rps( filebase,prf,pulseInCpi ,pulse_replica,ranges,track,fullarrayGeometry, 'windowrange',{@chebwin,40},'windowdoppler',{@chebwin,80},'nofft',true);
[sinr,dopplerbins,Rx]= theorysinr(trackt,fullarrayGeometry ,prf,pulseInCpi,lamda,reflextivityt,transmitterVelocity,receiverVelocity);
sv = exp(1i*2*pi*1/lamda * fullarrayGeometry*[cosd(angle);sind(angle);0]);
Vd = dftmtx(pulseInCpi );
pulse=((0:(pulseInCpi-1))')/prf;
gates = dopplerbins;
z = repmat(pulse,1,length(gates));
g = repmat(gates,length(pulse),1);
%%
doppler_sv = exp(1i * 2 * pi * z.*g);

BOM = kron(sv ,doppler_sv);

data = squeeze(Pn{1}(1,:,:,1));

%%
data = data(:);
w = repmat( data,1,size(BOM,2));
%%
sv = BOM;
sinrconv = abs(dot(w,sv)).^2./abs(dot(w,Rx*w));
sinrconv =  sinrconv/(pulseInCpi * size(fullarrayGeometry,1));
hold on;
plot(vc*dopplerbins,10*log10(abs(sinrconv )),'k','linewidth',3)
grid on
%%
Pn = Pn{1};
%for angle = [-90:10:90]
    angle = 0;
    sv = exp(-1i*2*pi*1/lamda * fullarrayGeometry*[cosd(angle);sind(angle);0]);
    for ii = 1:size(Pn,4)
        ii
        hold off;
        colormap jet;
        Px = Pn(:,:,:,ii);
        Px = reshape(Px,size(Px,1)*size(Px,2),size(Px,3))*sv;
        Px = reshape(Px,size(Pn,1),[]);
        zc = fftshift(20*log10(abs(sum(Px ,3).')),2);
        figure(1)
        imagesc(zc(200:400,:));
         caxis([-120 -40]);
        colorbar;
       
%         figure(2)
%         hold on;
%         plot(zc(290,:))
%         
%          drawnow;
%            figure(3)
%            hold on;
   % plot(zc(:,4))
    % caxis([20 90]);
    drawnow
    end