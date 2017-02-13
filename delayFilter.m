function [timeseries extra x] = delayFilter(transmittedTimeseries, state, clutterVoltage,...
    L, P, track ,Fs, Fc,C,...
    tStart,tEnd,antPos,parameters,first,Tz,rangeCorrect,sensorResponse,transmitGeometry,transmitOrientation,transmitShading,beamPower)
% DELAYFILTER FFT delay filter
% transmittedTimeseries -- function pointer to function that returns the
%                          timeseries
% state              -- state structure for timeseries function
% clutterVoltage    -- sqrt of clutter power (complex return for each
%                       clutter patch)
% L                  -- Length of fft in filter
% P                  -- Number of filter taps
% track              -- Tracks of transmitter, receiver and scatterers
%                        as a 4xMxN matrix the first four  rows are x,y,z
%                        and time.  M is the number of entries in the track
%                        N is the tracks N=1 is the track of the
%                        transmitter.  N=2 is the track of the receiver
%                        N=3 and greater is the track of each scatterer.
%                        The simulation will trace the track from the
%                        transmitter to the scatterer and then to the
%                        receiver
% Fs                 -- Sampling frequency in Hertz
% Fc                 -- Carrier frequency in Hertz
% C                  -- Propogation speed
% tStart             -- time of start of timeseries
% tEnd               -- time of end of timeseries
% antPos             -- Pos vector of anntenna location
% parameters         -- Parameters to pass to reference time series function
% first              -- Flag if first time code is called
% Tz                 -- Fixed sample offset (to decrease filter length)
% rangeCorrect      -- true if range attenuation is in effect
% sensorResponse    -- if not empty, contains a structure that defines the
%                       element response pattern
% transmitgeometry   -- sensorx3 geometry of transmit array
% transmitorientation -- orientation of transmitter  as (roll,pitch,yaw,t) matrix
% transmitShading    -- sensorx3 matrix of shading for transmit array
% select GPU or cpu processing by setting flag in code
CUDA = true;
ISFAST = false;
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
% set if real time series (non-baseband processing)
if(Fc==0)
    realTimeSeries = true;
else
    realTimeSeries = false;
end
% if(FARFIELD)
%     Tz = 0;
% end
% load time series
timess = [tStart:1/Fs:tEnd]-Tz;
state.create = false;
x = transmittedTimeseries(timess,state,parameters,Fs);

x = (x(:));

%if first call pad with zeros
if(first)
    x =[zeros(P/2-1,1);x];
end

xi = cast(single(x));
N = length(x);

D = L-P+1;

% calculate constants
if(realTimeSeries)
    waveNumber = cast(single(-1i * 2 * pi  *[0:1:L/2-1 0] * (Fs/L)).');
    myFwdFFT = @rfft;
    myInvFFT = @rifft;
    Kfc =0;
else
   waveNumber = cast(ifftshift(single(-1i * 2 * pi  *(-L/2:L/2-1) * (Fs/L)).'));
   %  waveNumber = cast(single(-1i * 2 * pi  *[0:1:L/2-1 0] * (Fs/L)).');
   % waveNumber = [waveNumber; waveNumber(end-1:-1:2,:)'.'];
    myFwdFFT = @fft;
    myInvFFT = @ifft;
    Kfc = -1j * 2 * pi * Fc;
end
c = 1;
totalLength = 0;

% break input into fft sections using overlap-save
endl = fix(N/D-1);
if(endl*D+L>length(xi))
    endl = endl-1;
end
Xrr = cast(casttype(zeros(L,endl+1)));
for r = 0:endl
    xr = xi(r*D + 1:r*D+L);
    Xrr(:,r+1) = xr;
    totalLength = totalLength + (L-fix(P/2)-1) - (fix(P/2)+c)+1;
    c = 0;
end
Xrr = myFwdFFT(circshift(Xrr,-length(Xrr)/2));


%sigma = cast(ones(size(sourcePos,2)));
if(first)
    t = [0:endl] * D/Fs +tStart ;
else
    t = [-1:endl] * D/Fs +tStart ;
end

%get delay ranges and power level of each source
if(isempty(beamPower{1}))
    [ ranges, clutterVoltage] =  ranginginformation(track,t+ (D/Fs)/2,cast,...
        rangeCorrect,sensorResponse,transmitGeometry,transmitOrientation,...
        transmitShading,clutterVoltage,Fc,C);
    %H = arrayfun(@delayFilterFreq,waveNumber,ranges.'/C,clutterVoltage ,Kfc,cast,Tz);
    %tic
    H = delayFilterFreq(waveNumber,ranges.'/C,clutterVoltage ,Kfc,cast,Tz);
    %toc
else
    [ ranges, clutterVoltage] =  ranginginformationtrim(track,t+ (D/Fs)/2,cast,...
        rangeCorrect,sensorResponse,transmitGeometry,transmitOrientation,...
        transmitShading,clutterVoltage,Fc,C,beamPower);
    H = cast(casttype(zeros(size(waveNumber,1),size(clutterVoltage,1))));
    ds = 20*log10(abs(clutterVoltage));
 
   
    for ii = 1:size(ds,1)
           minsnr = -100;
        if(max(ds(ii,:))<minsnr)
            minsnr =max(ds(ii,:))-20;
        end
        [jk I] = find(ds(ii,:)>minsnr);
        if(~isempty(I))
            rIter = ranges(I,ii);
            clutterVoltageIter = clutterVoltage(ii,I);
        else
            rIter = ranges(:,ii);
            clutterVoltageIter = clutterVoltage(ii,:);
        end
        H(:,ii) = delayFilterFreq2(waveNumber,cast(casttype(rIter.'/C)),cast(casttype(clutterVoltageIter)) ,casttype(Kfc),cast,casttype(Tz));
    end
    
end



H = myInvFFT(H);

win = fftshift(ones(P,1));

window = [bsxfun(@times,ones(fix(P/2),size(H,2)),win(1:fix(P/2)));...
    zeros(L-P+(mod(P,2)==1),size(H,2));...
    bsxfun(@times,ones(fix(P/2),size(H,2)),win(fix(P/2)+1:end))];
H = myFwdFFT(H.*cast(casttype(window) ));

index = fix(P/2):L-fix(P/2);


if(~ISFAST)
    if(first)
        Yj = Xrr .* H;
        yj =  myInvFFT( Yj);
        H2 = H(:,[1 1:end-1]);
        yj2 =  myInvFFT ( Xrr .*H2);
    else
        Yj = Xrr .* H(:,2:end);
        yj =  myInvFFT( Yj);
        H2 = H(:,[1:(end-1)]);
        yj2 =  myInvFFT ( Xrr .*H2);
    end
    
    
    
    blendLength = length(index);
    a = ([1:blendLength]'-1)/blendLength;
    y  = bsxfun( @times,a,yj(index,:)) + bsxfun( @times,(1-a),yj2(index,:));
else
    
    if(first)
        Yj = Xrr .* H;
    else
        Yj = Xrr .* H(:,2:end);
    end
    yj = myInvFFT( Yj);
    y  = yj(index,:);
end

timeseries = inverse_cast(gather(y));
timeseries = timeseries(2:end);

in = length(timeseries);
extra = N - in(1) -1;
timess = [tStart:1/Fs:tEnd];
state.create = false;

x = transmittedTimeseries(timess,state,parameters,Fs);
x = x(:);
if(first)
    x =[zeros(P/2-1,1);x];
end
x = x(P/2+1:end);
x = x(1:length(timeseries));

function [H] = delayFilterFreq2(K,ranges,sigma,Kfc,cast,Tz)
step_size =512;
%get index
start_index = 1;
end_index = step_size;
if(end_index>size(ranges,2))
    end_index = size(ranges,2);
end
index = start_index:end_index;
ranges = ranges-Tz;
% if(~(size(ranges,1)==1))
%     K = repmat(K,1,size(ranges,1));
% end

H = cast(zeros(size(K,1),size(ranges,1)));
while (start_index<=length(ranges))
    r = ranges(:,index);
    s = (sigma(:,index));
    M= (s.*exp(Kfc*r));
    M = repmat(M(:).',size(K,1),1);

    H = H + sum(reshape(exp( K* r(:).').*M,size(K,1),size(r,1),size(r,2)),3);
    
    % update index
    start_index = start_index + step_size;
    end_index = end_index + step_size;
    index = start_index:end_index;
    if(end_index>length(ranges))
        index = start_index:length(ranges);
    end
end

function [H] = delayFilterFreq3(K,ranges,sigma,Kfc,cast,Tz)
step_size =1000;
%get index
start_index = 1;
end_index = step_size;
if(end_index>size(ranges,2))
    end_index = size(ranges,2);
end
index = start_index:end_index;
% if(~(size(ranges,1)==1))
%     K = repmat(K,1,size(ranges,1));
% end

H = cast(zeros(size(K,1),size(ranges,1)));

while (start_index<=length(ranges))
    s = (sigma(:,index));
    sig = single(s);
    r = ranges(:,index);
    r = r-Tz;
    %if(size(ranges,1)==1)
     %   H = exp( K* r)*M;
   % else
        for ii = 1:size(ranges,1)
            I = find(abs(sig(ii,:))>.01*max(abs(sig(ii,:))));
            s2 = s(ii,I);
            r2  = r(ii,I);
            M= (s2.*exp(Kfc*r2)).';
            H(:,ii) = H (:,ii)+ exp( K* r2)*M;
        end
   % end
    % update index
    start_index = start_index + step_size;
    end_index = end_index + step_size;
    index = start_index:end_index;
    if(end_index>length(ranges))
        index = start_index:length(ranges);
    end
end