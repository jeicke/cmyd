% Test of proper scattering range
% tests that we have the poper range to scatterer
% RADAR PARAMETERS
NUMBEROFSCATTERPOINTS = 20;
Fs = 60E6; %samples/second
Fc = 9410E6; %hertz
rpm = 24; %rotations per minute of radar
verticalBeamwidth = 25; % in degrees
horizontalBeamwidth = 1.8; % in degrees
prf = 3000; % in hz
pulseBandwidth = 30E6; % in hz
pulseLength = 2400E-9; % in seconds
C = 299792458; %m/s
filebase = 'C:\Users\rnl\Desktop\testdata\test';
image_filename = 'C:\Users\rnl\Desktop\ISAR code\test3.png';
invert_image = true;
mask;
% source and receiver locations
receiverPosition = [0 0 0];   
transmitterPosition = [2000 1000 40];
receiverVelocity = [30 0 0];   
transmitterVelocity= [0 0 0];

% derived parameters
lamda = C/Fc;
arrayHeight = (1.27 * lamda)/(verticalBeamwidth *pi/180);
arrayLength = (1.27 * lamda)/(horizontalBeamwidth *pi/180);
dwell_time = 1/(rpm/60) * horizontalBeamwidth/360;
Npulses = floor(dwell_time * prf);%number of pulse in look direction
plength = round(Fs/ prf);
% just do a few look for this simulation
time = 10.1/(rpm/60);

% set track for rotating transmit antenna
t = [0:dwell_time:time];
orientation = zeros(4,length(t));
orientation(3,:) = [0:horizontalBeamwidth:(length(t)-1)*horizontalBeamwidth],360;% - radar_velocity(3)/2  * time;
orientation(4,:) = t;

%geometry of transmit array
[geometry w h]= loadgeometry('rectangle',[],arrayLength,arrayHeight,lamda);

% make transmitter and receiver tracks
[ track] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );

% make scattered point pegged to transmitter (like we're on a boat)

clear reflextivity;

clear reflextivity;
for ii = 1:size(scatteredpoint,2);
    %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    
    reflextivity(ii) = sigma(ii);
    track(:,:,ii+2) = track(:,:,1);
    track(1,:,ii+2) =track(1,:,ii+2) +scatteredpoint(1,ii);
    track(2,:,ii+2) =track(2,:,ii+2) +scatteredpoint(2,ii);
    track(3,:,ii+2) =track(3,:,ii+2) +scatteredpoint(3,ii);
end
w =chebwin(40,70);
w = kron(w,[1 1 1]');
w = w/sum(abs(w));
%%
radsim(filebase,time,track,0,...
       'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
       'fs',Fs,'fc',Fc,'transmitgeometry',geometry,...
       'transmitorientation',orientation,'clutter',reflextivity',...
       'transmitshading',w); 
   %%
pulse_replica=m_chirp(pulseBandwidth,pulseLength,0,Fs)'.';
pulse_replica=pulse_replica(end:-1:1);
rngspa= 299792458/Fs;
ranges = [0:1:plength] * rngspa- (length(pulse_replica)-2) * rngspa;

rm = norm(receiverPosition -transmitterPosition);
rm = floor(rm/rngspa)+length(pulse_replica)-1
step = floor(1500 / rngspa);
rangeGates = [max([1 rm-step]):1:min([rm+step length(ranges)])];

[D blockTimes] = processradar( filebase,prf,1,pulse_replica,rangeGates);
ranges = ranges(rangeGates );
dt = 1/prf;
t = ([0:1:(size(D,3)-1)])*(dt * size(D,2))+dt/2;
%%
pixelSpaceing = 2;
x = 1:200;
x = (x-mean(x)) * pixelSpaceing ;

y = 1:200;
y = (y-mean(y)) * pixelSpaceing ;

z = ones(size(y)) * -40;
clear trialpoint;
count = 1;
for ii = 1:length(x)
    for jj = 1:length(y)
        
        
        trialpoint(:,count) = [x(ii) y(jj) z(jj)];
        count = count + 1;
    end
end
[ trackt] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );
for ii = 1:size(trialpoint,2);
    %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    
    
    trackt(:,:,ii+2) = trackt(:,:,1);
    trackt(1,:,ii+2) =trackt(1,:,ii+2) +trialpoint(1,ii);
    trackt(2,:,ii+2) =trackt(2,:,ii+2) +trialpoint(2,ii);
    trackt(3,:,ii+2) =trackt(3,:,ii+2) +trialpoint(3,ii);
end
for jj = 1:100
    jj
    az = interp1(orientation(4,:),orientation(3,:),blockTimes(jj));

    
    [r_transmitter_scatter aztransmit detransmit]= computeRange(trackt, blockTimes(jj),1);
    
    % next comput range from scatterer to receiver and the angles
    %bsxfun(@plus,scatterPos,-antPos)
    [r_scatter_ant azreceive dereceive]= computeRange(trackt, blockTimes(jj),2);
    rangelook = (zeros(size(r_transmitter_scatter)));
    
    rangelook(1,:) = (r_transmitter_scatter(1,:));
    
    % add path from transmitter to scatter and from scatter to receiver  together for scattered path
    if(size(r_transmitter_scatter,1)>1)
        rangelook(2:end,:) = (r_transmitter_scatter(2:end,:) + r_scatter_ant(2:end,:));
    end
    aztransmit = aztransmit * 180/pi;
    aztransmit(aztransmit<0 ) = aztransmit(aztransmit<0 ) + 360;
    aztransmit = mod(aztransmit,360);
    
    
    
    s = interp2(az,ranges,abs(p),aztransmit(2:end),rangelook(2:end),'cubic');
    s = reshape(s,200,[]);
end
%%
figure;
M = double(10*log10(abs(s)));
[xx yy] = meshgrid(x,y);
pcolor(xx,yy,M)
xlabel('POSITION (m)')
ylabel('POSITION (m)')
caxis([max(max(M))-30 max(max(M))])
shading flat
axis equal
% hold on;plot(scatteredpoint(1,:),scatteredpoint(2,:),'ko','MarkerSize',16, 'MarkerFaceColor','k')
figure;
pcolor(xx,yy,reshape( 15 * floor(rangelook(2:end)/30),200,[]))
xlabel('POSITION (m)')
ylabel('POSITION (m)')
shading flat
axis equal

figure;plottrack(track)

