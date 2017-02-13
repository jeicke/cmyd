%% now staggered pri postDoppler
snapshots = 20;
guardband = 1;
sigma = 1/100;
[Pn] =rps( filebase,prf,pulseInCpi ,pulse_replica,ranges,track,fullarrayGeometry, 'windowrange',{@chebwin,40},'windowdoppler',{@chebwin,80});
[Pf] =rps( filebase,prf,pulseInCpi ,pulse_replica,ranges,track,fullarrayGeometry, 'windowrange',...
    {@chebwin,40},'windowdoppler',{@chebwin,80},'nofft',true);
Pn = Pn{1};
Pf = Pf{1};
clear P;
wd = chebwin(pulseInCpi,80);
wd = wd./norm(wd);
for angle = [5]
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
            for dpbin = 2%1:size(Px,1)
                
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