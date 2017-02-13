[ trackt] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );
for ii = 1:size(trialpoint,2);
    %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    trackt(:,:,ii+2) = trackt(:,:,1);
    trackt(1,:,ii+2) = trackt(1,:,ii+2)+ trialpoint(1,ii);
    trackt(2,:,ii+2) = trackt(2,:,ii+2)+ trialpoint(2,ii);
    trackt(3,:,ii+2) = trackt(3,:,ii+2) +trialpoint(3,ii);
end
sourceVelocity = transmitterVelocity;
for ii = 1:(size(trackt,3)-2)
 p1 = trackt(1:3,1,1)-trackt(1:3,1,ii+2);
 p2 = trackt(1:3,1,ii+2)-trackt(1:3,1,2);
 t1 = norm(p1);
 t2 = norm(p2);
 u1 = p1/t1;
 u2 = p2/t2;
 
 trrange = t1+t2;
 %p = [x2;y(ii,1);z(ii,1)] /norm( [x2;y(ii,1);z(ii,1)]);
 d1=  (sourceVelocity-transmitterVelocity) * u1* (Fc)/299792458;
 d2=  (receiverVelocity-sourceVelocity) * u2* ( Fc)/299792458;
 doppler(ii) = d1 + d2;
 %plot([doppler doppler],[ranges(1) ranges(end)]/1000,'r--');
end
%%
%[D2 ] = pulsecompression(D2,pulse_replica,'cheb',30);

P = zeros(1,length(x)*length(y),1);

%for ii = 1:size(D2,3)
ii = 1;
az = interp1(orientation(4,:),orientation(3,:),blockTimes(ii:ii+size(D2,3)));
az = mod(azx(1+(ii-1) * frameLength:ii*frameLength),360);




[r_transmitter_scatter aztransmit detransmit pos1]= computeposition(trackt,  blockTimes(1)',1);
aztransmit = aztransmit * 180/pi;
aztransmit(aztransmit<0 ) = aztransmit(aztransmit<0 ) + 360;
aztransmit = mod(aztransmit,360);
% next comput range from scatterer to receiver and the angles
%bsxfun(@plus,scatterPos,-antPos)
%[r_scatter_ant1 azreceive dereceive pos1]= computeposition(trackt,  blockTimes(1)',2);
r_scatter_ant = zeros(size(r_transmitter_scatter));
for jj  = 1:size(trialpoint,2)
    [h ind] = min(abs( aztransmit(jj)-az));
    [temp j1 j2 pos2]= computeposition(trackt(:,:,[1 2 jj+2]),  blockTimes(ind)',2);
    r_scatter_ant(jj+1) = temp(2);
    
    
    p1 = pos1(:,jj+2);
    p2 = pos2(:,3);
    t1 = norm(p1);
    t2 = norm(p2);
    u1 = p1/t1;
    u2 = p2/t2;
    
    trrange = t1+t2;
    %p = [x2;y(ii,1);z(ii,1)] /norm( [x2;y(ii,1);z(ii,1)]);
    d1=  (sourceVelocity-transmitterVelocity) * u1* (Fc)/299792458;
    d2=  (receiverVelocity-sourceVelocity) * u2* ( Fc)/299792458;
    doppler(jj) = d1 + d2;
end
rangelook = (zeros(size(r_transmitter_scatter)));

rangelook(1,:) = (r_transmitter_scatter(1,:));

% add path from transmitter to scatter and from scatter to receiver  together for scattered path
if(size(r_transmitter_scatter,1)>1)
    
    rangelook(2:end,:) = (r_transmitter_scatter(2:end,:) + r_scatter_ant(2:end,:));
end



M = fftshift(abs(D2),2);

s = interp3(f,ranges,az,M,doppler(1:end)',rangelook(2:end),aztransmit(2:end),'linear',eps);

P2(1,:,ii) =s;

%end
