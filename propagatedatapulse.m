function L = propagatedatapulse(transmittedTimeseries, track, clutterPower,Fs,Fc,C,run_length,...
    fileNameBase,maxSystemDelay,minSystemDelay,...
    parameters,geometry,rangeCorrect,sensorResponse,transmitGeometry,transmitOrientation,transmitShading)
%sample rate


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
%
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

while iBlocks<totalBlocks
    % delay time series
   
    tStart = index/Fs+Tz;
    tEnd = tStart+L/Fs;
    
    if(tEnd>run_length)
        tEnd =run_length;
    end

    times = [tStart:1/Fs:tEnd];


    state.maxSystemDelay = maxSystemDelay+2*L/Fs;
    state.create = true;
    [~, state] = transmittedTimeseries(times(extraFilter+2:end),state,parameters,Fs);
    clear TimeSeries;
    % calculate received time series from all sources
    for ii = 1:size(geometry,1)
         pos(:,1) = geometry(ii,:)';
         if(exist('sensor_response','var')&&~isempty(sensorResponse))
             sensorResponse.uc =  sensorResponse.alluc(ii,:);
         end
        [TimeSeries(ii,:),  extraFilter, data] = delayfilterblock(transmittedTimeseries, state, clutterPower,L,...
            track,Fs,Fc,C,...
            tStart,tEnd,geometry(ii,:)',parameters,...
            firstFlag,rangeCorrect,sensorResponse,transmitGeometry,transmitOrientation,transmitShading,beamPower,Tz);
        
        
        %add noise
        %totalTimeSeries= addNoise(totalTimeSeries,power_per_sample,sim);
        
        % write delayed timeseries
    
        
    end
    TimeSeries =  TimeSeries(:,1:L- floor(length(pulse)/2));
     fp1 = TestWrite(TimeSeries(:).',firstFlag,sprintf('%s_%s',fileNameBase,'received.bin'),Fs,Fc,size(geometry,1),fp1);
     if(firstFlag)
         fp2 = TestWrite(double(data).',firstFlag,sprintf('%s_%s',fileNameBase,'transmitted.bin'),Fs,Fc,1,fp2);
     end
     % write transmitted timeseries
    firstFlag = false; 
    
    
    iBlocks = iBlocks + 1;
    index = index +  blockLength;
    t = toc;
    if( (t >5)||iBlocks ==1)
        
        fprintf(1,'\tPercent done: %2.2f\n',min([100 100*(iBlocks)/totalBlocks]));
        tic
    end
    
end
L = L- floor(length(pulse)/2);
fprintf(1,'\tPercent done: %2.2f\n',min([100 100*(iBlocks)/totalBlocks]));
fclose(fp1);
fclose(fp2);
