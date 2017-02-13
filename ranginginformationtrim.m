function [ ranges clutterVoltage] = ranginginformation(track,times,cast,...
    rangeCorrect,sensorResponse,transmitGeometry,transmitOrientation,...
    transmitShading,clutterVoltage,Fc,C,beamPower)
%COMPUTERANGES computes all ranges for track points at given times
%   Detailed explanation goes here
% add antenna position to track

az = [0:.1:360-.1];
de = [-90:.1:89];

% for scatter the path is from transmitter to scatterer to receiver
% first compute path from transmitter to scatterer and the angles to
% the scatterer
FARFIELD = false;
[r_transmitter_scatter aztransmit detransmit]= computeRange(track, times,1);

% next comput range from scatterer to receiver and the angles
%bsxfun(@plus,scatterPos,-antPos)
[r_scatter_ant azreceive dereceive]= computeRange(track, times,2);
ranges = (zeros(size(r_transmitter_scatter)));
rangeAttenuationTransmitterScatter = ranges;
rangeAttenuationScatterReceiver = ranges;
rangeAttenuationTransmitterScatter(2:end,:) = 1./r_transmitter_scatter(2:end,:);
rangeAttenuationScatterReceiver(2:end,:) = 1./r_scatter_ant(2:end,:);

% direct path from transmitter to receiver
ranges(1,:) = (r_transmitter_scatter(1,:));
rangeAttenuationTransmitterScatter(1,:) = 1./r_transmitter_scatter(1,:);
rangeAttenuationScatterReceiver(1,:) = 1;
% add path from transmitter to scatter and from scatter to receiver  together for scattered path
if(size(r_transmitter_scatter,1)>1)
    ranges(2:end,:) = (r_transmitter_scatter(2:end,:) + r_scatter_ant(2:end,:));
    rangeAttenuationTransmitterScatter(2:end,:) = 1./r_transmitter_scatter(2:end,:);
    rangeAttenuationScatterReceiver(2:end,:) = 1./r_scatter_ant(2:end,:);
end

if(~isempty(transmitGeometry))
    
     azt = mod(interp1q(transmitOrientation(4,:)',transmitOrientation(3,:)', (times)'),360);
 
%     A1 = sin(aztransmit);
%      A2 = cos(aztransmit);
%      B1 = sind(azt);
%      B2 = cosd(azt);
    for ii = 1:size(aztransmit,2)
       % rot = [cosd(azt(ii))  -sind(azt(ii)) 0; sind(azt(ii))  cosd(azt(ii)) 0; 0 0 1];
        if(~iscell(transmitGeometry))
            geometry = transmitGeometry * rot;
            [power] = beampattern(Fc,180/pi * aztransmit(:,ii),...
                180/pi*detransmit(:,ii),ones(size(geometry,1),1),[], transmitShading,geometry,[],C);
        else
           % geometry = transmitGeometry{1} * rot;
            %[powerh] = beampattern(Fc,180/pi * aztransmit(:,ii),...
         %       zeros(size(aztransmit(:,ii))),ones(size(geometry,1),1),[], transmitShading{1},geometry,[],C);
            azs = 180/pi *  aztransmit(:,ii);
            azs(azs<0) = azs(azs<0) + 360;
            azs = mod(azs,360);
            azi = az+azt(ii);
             azi(azi<0) = azi(azi<0) + 360;
             azi = mod(azi,360);
            powerh = interp1(azi',sqrt(beamPower{1}),azs);
            geomv = transmitGeometry{2};
            powerv = interp1q(de' ,sqrt(beamPower{2}),180/pi*detransmit(:,ii));
%             [powerv] = beampattern(Fc,180/pi * aztransmit(:,ii),...
%                 180/pi*detransmit(:,ii),ones(size(geomv ,1),1),[],transmitShading{2} ,geomv ,[],C);
            power = abs(powerh.*powerv);
        end
         beamAttenuation(:,ii)  = power;
%          s = [A1(:,ii) A2(:,ii)];
%          m = [B1(ii) B2(ii)];
%          m2 = repmat(m',1,size(s,1));
%          DIFF = acos( dot(s.',m2));
%          for jj = 1:size(aztransmit,1)
%              if(aztransmit(jj,ii)<0)
%                  aztransmit(jj,ii)  = aztransmit(jj,ii)  + 2*pi;
%              end
%             
% 
%            %  DIFF = acos(dot([A1(jj,ii) A2(jj,ii)],[B1(ii) B2(ii)]));
%               %fprintf(1,'TRANSMITE BEAM %3.1f ANGLE TO SCATTER %3.1f DIFF %3.1f FLAG %d\n',azt(ii),aztransmit(jj,ii) * 180/pi,DIFF * 180/pi,DIFF>pi/2)
%              if(DIFF(jj)>pi/2)
%                 
%                  beamAttenuation(jj,ii) = beamAttenuation(jj,ii)*10^(-40/20);
%              end
%          end
    end
    
else
    
end
 beamAttenuation = ( beamAttenuation);
clutterVoltage = (clutterVoltage);
rangeAttenuationTransmitterScatter = (rangeAttenuationTransmitterScatter);
rangeAttenuationScatterReceiver = (rangeAttenuationScatterReceiver);
clutterVoltage = repmat(clutterVoltage,size(ranges,2),1);
if(rangeCorrect)
    clutterVoltage = clutterVoltage.*rangeAttenuationTransmitterScatter.';
    clutterVoltage = clutterVoltage.*rangeAttenuationScatterReceiver.';
end
if(~isempty(transmitGeometry))
    clutterVoltage = clutterVoltage .*  beamAttenuation';
end


% if we want to do a farfield computation...
if(FARFIELD)
    ranges = scatterPos' * [sin(azreceive).*cos(dereceive);cos(azreceive).*cos(dereceive);sin(dereceive)];
end
if(exist('sensorResponse','var')&&~isempty(sensorResponse))
    crfun = sensorResponse.bpfun;
    Cr = (crfun(azreceive, dereceive, sensorResponse)');
    clutterVoltage = clutterVoltage.*Cr;
end

end

