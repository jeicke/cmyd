[ trackt] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );
for ii = 1:size(trialpoint,2);
    %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    trackt(:,:,ii+2) = trackt(:,:,1);
    trackt(1,:,ii+2) = trackt(1,:,ii+2) + trialpoint(1,ii);
    trackt(2,:,ii+2) = trackt(2,:,ii+2) + trialpoint(2,ii);
    trackt(3,:,ii+2) = trackt(3,:,ii+2) + trialpoint(3,ii);
end
%[D2 ] = pulsecompression(D2,pulse_replica,'cheb',30);

P = zeros(1,length(x)*length(y),size(D2,3));

for ii = 1:size(D2,3)
    ii

        az = interp1(orientation(4,:),orientation(3,:),blockTimes(ii:ii+size(D2,2)));
        az = mod(azx(1+(ii-1) * frameLength:ii*frameLength),360);
        
        
        
        
        [r_transmitter_scatter, aztransmit, detransmit]= computeRange(trackt,  blockTimes2(1,ii)',1);
        aztransmit = aztransmit * 180/pi;
        aztransmit(aztransmit<0 ) = aztransmit(aztransmit<0 ) + 360;
        aztransmit = mod(aztransmit,360);
        % next comput range from scatterer to receiver and the angles
        %bsxfun(@plus,scatterPos,-antPos)
        [r_scatter_ant1, azreceive, dereceive]= computeRange(trackt,  blockTimes2(1,ii)',2);
        r_scatter_ant = zeros(size(r_scatter_ant1));
        for jj  = 1:size(trialpoint,2)
             [h ind] = min(abs( aztransmit(jj)-az));
             [temp]= computeRange(trackt(:,:,[1 2 jj+2]),  blockTimes2(ind,ii)',2);
             r_scatter_ant(jj+1) = temp(2);
        end
        rangelook = (zeros(size(r_transmitter_scatter)));
        
        rangelook(1,:) = (r_transmitter_scatter(1,:));
        
        % add path from transmitter to scatter and from scatter to receiver  together for scattered path
        if(size(r_transmitter_scatter,1)>1)
           
            rangelook(2:end,:) = (r_transmitter_scatter(2:end,:) + r_scatter_ant(2:end,:));
        end
        
        
        
        M = abs(D2(:,:,ii));
        M(:,az>70&az<110) = median(M(:));
        s = interp2(az,ranges,M,double(aztransmit(2:end)),rangelook(2:end),'linear',eps);
        P(1,:,ii) =s;

end
