for mm = 1
snapshots = 30;
guardband = 1;
if(mm == 1)
    igma = .0001;
    pulseInCpi = 32;
    staggers =1;
    label = 'conv';
end
if(mm == 2)
    sigma = .001;
    pulseInCpi = 64;
    staggers =1;
    label = 'ABF';
end
if(mm == 3)
    sigma = .0001;
    pulseInCpi = 64;
    staggers =2;
    label = 'STAGGER';
end
sigma1 = .001;


[Pn] =rps( filebase,prf,pulseInCpi ,pulse_replica,ranges,track,fullarrayGeometry, 'windowrange',{@chebwin,40},'windowdoppler',{@chebwin,80},'postdoppler',staggers);
[Pf] =rps( filebase,prf,pulseInCpi ,pulse_replica,ranges,track,fullarrayGeometry, 'windowrange',...
    {@chebwin,40},'windowdoppler',{@chebwin,80},'nofft',true);
% make temporal steering vector matrix
N = size(Pn{1},1);
dsv = dftmtx(N);
wd =chebwin(pulseInCpi,80);
clear P Px;
clear sinr
for angle = [90]
    %spatial steering vector
    sv = exp(1i*2*pi*1/lamda * fullarrayGeometry*[sind(angle);cosd(angle);0]);
    Vd = dftmtx(pulseInCpi );
    BOM = kron(sv ,Vd);
    
    In = eye(size(fullarrayGeometry,1));
    clear P P2;
    for ii = 1%:size(Pn{1},4)
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
        Pxf = Pf{1}(:,:,:,ii);
        count = 1;
        for rg = 200:400
            
            for dpbin = 1:(size(Px{1},1)-staggers +1)
                
                BOM = kron(sv ,Vd(:,dpbin ));
                % estimate spatial-doppler correlation matrix for each range using N snapshots
                % around given range value with a given guard band
                startLeft = max([rg - (snapshots+guardband) 1 ]);
                endLeft = max([rg  - guardband-1 1]);
                startRight = min([rg  + guardband+1 size(Px{1},2)-1]);
                endRight = min([rg +(snapshots+guardband) size(Px{1},2)-1]);
                s = [];
                data =Pxf(rg,:);
                Fm = zeros(length(dsv),numel(Pn));
                
                for pdb =  1:numel(Pn)
                    if(numel(sv)~=1)
                        s = [s squeeze((Px{pdb}(dpbin,[startLeft:endLeft startRight:endRight],:)))'.'];
                    else
                         s = [s squeeze((Px{pdb}(dpbin,[startLeft:endLeft startRight:endRight],:)))'];
                    end
                    endIndex = pulseInCpi -(numel(Pn)-pdb);
                    Fm(pdb:endIndex,pdb) = dsv(:,dpbin);
                end
                T = kron(In,Fm.*wd);
                
                R = (s' * s)/size(s,1);
                R = R + mean(diag(R))*sigma*eye(size(R));
                if(mm == 1)
                    w = (T'*BOM);

                else
                    w = inv(R)*(T'*BOM);
                    
                    w = w./dot(w,R*w);
                end
                v = [];
                for pdb =  1:numel(Pn)
                    v = [v squeeze((Px{pdb}(dpbin,rg,:)))'];
                end
                v = v.';
                wa = T*w;
                vall = BOM;
                P(count,dpbin) = abs(w'*v);
                Ri = data'*data;
                Ri = Ri+mean(diag(Ri))*sigma1*eye(size(Ri));
                wa = wa.';
                P2(count,dpbin) =abs(wa*Ri*wa');
                
                sinr(count,dpbin) = abs(wa'.'*vall).^2./abs(wa*Ri*wa').^2;
            end
            count = count + 1;
        end
%         figure(8)
%         plot(20*log10(abs(Vd*wa)))
%         hold on
%         plot(20*log10(abs(Vd*data.')));
%         figure(4)
%         colormap jet
%         imagesc(doppler_bins,[1 12.5],12.5+fftshift(20*log10(abs(P)),2))
%         colorbar
%         %caxis([-120 -40]);
%         drawnow;
        figure
        colormap jet
        p = fftshift(10*log10(abs(P2)),2);
        imagesc(doppler_bins,[1 12.5],p)
        colorbar
        caxis([max(p(:))-80 max(p(:))]);
        drawnow;
        title(label);
    end
end
sinrAll{mm} = 20*log10(fftshift(mean(sinr(:,:),1)));
end
%%
doppler_bins =vc*linspace(-prf/2,prf/2,pulseInCpi);
figure;
plot(doppler_bins,sinrAll{1}-sinrAll{1}(1))
hold on;
plot(doppler_bins,sinrAll{2}-sinrAll{2}(1))
plot(doppler_bins(1:end-2),sinrAll{3}-sinrAll{3}(1))
grid on;
legend('CONVENTIONAL','ADAPTIVE BF','STAGGER STAP')
%%
[r_scatter_ant, azreceive, dereceive]= computeRange(track, 0*pulseInCpi *1/prf,2);
figure
az = -90:.1:90;
el = zeros(size(az));
for ii = 280
    x = squeeze(Pn{1}(2,ii,:,1));
    
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