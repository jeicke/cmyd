resolution = 1;
% AZ = ( ([-90:resolution:90]) * pi/180);
 %AZ = ( ([0:resolution:180]) * pi/180);
 % AZ = ( ([40:resolution:61]) * pi/180);
  AZ = ( (90+[-45:45]) * pi/180);
rngspa= 299792458/(Fs);

clutter_range = [1 12.5 ];
ranges = 1+[clutter_range(1) * 1000:rngspa:(clutter_range(2))*1000];
%ranges = ranges(10:11);
% compute clutter area
[clutterarea ] = bistaticpatcharea(AZ,ranges,receiverPosition(3)/1000,gamma);
%ranges = ranges(13);
%figure(1)
%plot(ranges,clutterarea.*ranges*pi);
% make clutter patches
trackCenter = [0 0 0];
fprintf(1,'RANGE RESOLUTION %3.2f m\n',C/(2*pulseBandwidth))
fprintf(1,'CROSS-RANGE RESOLUTION %3.2f m\n',norm(trackCenter-receiverPosition)*lamda/(2*norm(receiverVelocity)*time))
clear scatteredpoint;
clear reflextivity;
count = 1;
% scatteredpoint(1,count) = 5000;
% scatteredpoint(2,count) = 0;
% scatteredpoint(3,count) = 0;
% reflextivity(count)  = clutterarea(ri);
for azi = 1:length(AZ)
    for ri = 1:length(ranges)
        scatteredpoint(1,count) = ranges(ri)*sin(AZ(azi));
        scatteredpoint(2,count) = ranges(ri)*cos(AZ(azi));
        scatteredpoint(3,count) = 0;
        reflextivity(count)  = clutterarea(ri);
        count = count + 1;
    end
end

for ii = 1:size(scatteredpoint,2)
    %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    track(:,:,ii+2) = track(:,:,1);
    track(1,:,ii+2) =trackCenter(1);
    track(2,:,ii+2) =trackCenter(2);
    track(3,:,ii+2) =trackCenter(3);
    %track(1:3,:,ii+2) = 0;
    track(1,:,ii+2) =track(1,:,ii+2)+scatteredpoint(1,ii);
    track(2,:,ii+2) =track(2,:,ii+2)+scatteredpoint(2,ii);
    track(3,:,ii+2) =track(3,:,ii+2)+scatteredpoint(3,ii);
    
end

%now just use test ranges
testRange = 1;
count = 1;
clear scatteredpointt;
clear reflextivityt;
clear trackt;
for azi = 1:length(AZ)
    for ri = testRange
        scatteredpointt(1,count) = ranges(ri)*sin(AZ(azi));
        scatteredpointt(2,count) = ranges(ri)*cos(AZ(azi));
        scatteredpointt(3,count) = 0;
        reflextivityt(count)  = clutterarea(ri);
        count = count + 1;
    end
end
trackt(:,:,1) = track(:,:,1);
trackt(:,:,2) = track(:,:,2);
for ii = 1:size(scatteredpointt,2)
    %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    trackt(:,:,ii+2) = track(:,:,1);
    trackt(1,:,ii+2) =trackCenter(1);
    trackt(2,:,ii+2) =trackCenter(2);
    trackt(3,:,ii+2) =trackCenter(3);
    %track(1:3,:,ii+2) = 0;
    trackt(1,:,ii+2) =trackt(1,:,ii+2)+scatteredpointt(1,ii);
    trackt(2,:,ii+2) =trackt(2,:,ii+2)+scatteredpointt(2,ii);
    trackt(3,:,ii+2) =trackt(3,:,ii+2)+scatteredpointt(3,ii);
    
end