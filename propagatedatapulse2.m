function L = propagatedatapulse(transmittedTimeseries, track, clutterPower,Fs,Fc,C,run_length,...
    fileNameBase,maxSystemDelay,minSystemDelay,...
    parameters,geometry,rangeCorrect,sensorResponse,transmitGeometry,transmitOrientation,transmitShading)
%sample rate

CUDA = true;

if(CUDA)
    %cast = @GPUsingle;
    cast = @gpuArray;
    casttype = @single;
    inverse_cast = @gather;
else
    cast = @single;
    casttype = @single;
    inverse_cast = @double;
end
firstFlag = true;

iBlocks = 0;
index = 0;
state = [];
%fixed delay
Z = floor((minSystemDelay)*Fs);
Tz = Z/Fs;

state.create = true;
state.maxSystemDelay = maxSystemDelay;
[~, state] = transmittedTimeseries(0,state,parameters,Fs);

state.create = true;
pulse = state.chirp(abs(state.chirp)>0);
L = floor((maxSystemDelay-(minSystemDelay))*Fs);
if(L<2*length(pulse ))
    L = 2*length(pulse );
end
L = 2.^(nextpow2(L));
%make transmitted pulse
state.create = false;
x = pulse;
x(abs(x)==0) = [];
x = (x(:));
xi = cast(casttype(x));
if(CUDA)
    x2 = gpuArray.zeros(L,1);
else
    x2 = zeros(L,1);
end
x2(1:length(xi)) = xi;
%x2 = circshift(x2,-length(xi)/2);
% x2 = zeros(size(x2));
% x2(1) = 1;
Xrr = fft(x2,L);

fprintf(1,'***FFT Length %d\n',L);
prf = parameters(3);
numberOfPulses = floor(run_length*prf);
blockLength = round(Fs*1/prf);
totalBlocks = numberOfPulses;

extraFilter = -1;

fp1 = [];
fp2 = [];
beamPower{1} = [];
if(1)
    % calculate beampatterns
    if(iscell(transmitGeometry))
        %calculate horizontal pattern
        az = [0:.1:360-.1];
        [beamPower{1}] = beampattern(Fc,az,...
            zeros(size(az)),ones(size(transmitGeometry{1},1),1),[], transmitShading{1},transmitGeometry{1},[],C);
        I = find(az>90&az<270);
        beamPower{1}(I) = beamPower{1}(I) *10^(-60/10);
        %calculate vertical pattern
        de = [-90:.1:89];
        [beamPower{2}] = beampattern(Fc,zeros(size(de)),...
            de,ones(size(geometry,1),1),[], transmitShading{2},transmitGeometry{2},[],C);
        
    end
end
tic
cpiTimeStep = round(Fs*1/prf)/Fs;
blockTimes =  [0:(numberOfPulses-1)]*cpiTimeStep+(maxSystemDelay-minSystemDelay)/2 +minSystemDelay;
blockStep =8;
iBlocks = 1;
while iBlocks<totalBlocks
    index = iBlocks:(iBlocks+blockStep);
    index(index>totalBlocks) = [];
    stepTimes = blockTimes(index);
    allTimeSeries= [];
    % calculate received time series from all sources
    
    for ii = 1:size(geometry,1)
        thisTrack = track;
        % shift receiver track to antenna position
        thisTrack(1,:,2) = thisTrack(1,:,2)+geometry(ii,1);
        thisTrack(2,:,2) = thisTrack(2,:,2)+geometry(ii,2);
        thisTrack(3,:,2) = thisTrack(3,:,2)+geometry(ii,3);
        [TimeSeries, data] = delayfilterblock2(Xrr , state, clutterPower,L,...
            thisTrack,Fs,Fc,C,stepTimes,[],parameters,...
            firstFlag,rangeCorrect,sensorResponse,transmitGeometry,transmitOrientation,transmitShading,beamPower,Tz);
        
        allTimeSeries(:,:,ii) = TimeSeries;
    end
    allTimeSeries = permute(allTimeSeries,[1 3 2]);
    fp1 = TestWrite(allTimeSeries(:).',firstFlag,sprintf('%s_%s',fileNameBase,'received.bin'),Fs,Fc,size(geometry,1),fp1);
    if(firstFlag)
        fp2 = TestWrite(double(data).',firstFlag,sprintf('%s_%s',fileNameBase,'transmitted.bin'),Fs,Fc,1,fp2);
    end
    % write transmitted timeseries
    firstFlag = false;
    iBlocks = iBlocks + blockStep+1;
    
    t = toc;
    if( (t >5)||iBlocks ==1)
        
        fprintf(1,'\tPercent done: %2.2f\n',min([100 100*(iBlocks)/totalBlocks]));
        tic
    end
    
end

fprintf(1,'\tPercent done: %2.2f\n',min([100 100*(iBlocks)/totalBlocks]));
fclose(fp1);
fclose(fp2);
