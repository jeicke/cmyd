function [P] =bpm( fileNameBase,prf,pulseReplica,ranges,track,geometry,xi,yi,zi,varargin)
% BPM Performs SAR time domain correlation
% Stands for:
% Back Propogation Method (BPM)
% perfoms backpropogation of SAR data
% INPUT
% 'fileNameBase'      -- base filename for data to process=
% 'prf'               -- Pulse repetition frequency in Hz
% 'pulseReplica'      -- IQ samples of transmitted fast time pulse
% 'track'             -- Track of simulation objects
% 'xi'                -- xi Evaluation sample points
% 'yi'                -- yi Evaluation sample points
% 'zi'                -- zi Evaluation sample points
% Keyword inputs
% 'cuda'              -- 'cuda','true' sets GPU computation to on
%                        default off
% 'windowrange'       -- 'windowrange','@hann' will use a hann range
%                         window, for example.  Default uniform
% 'windowcrossrange'  -- same thing but in the cross range diminision
%                     -- Default uniform
% 'timewindow'        -- length of time segment in seconds this is to
%                        implement Welch's method
%                        in the SAR world these techniques are known as
%                        'multiple-look processing'
%                        default set to one segment for entire image
% 'timeoverlap'       -- length of time in seconds that segments overlap
%                        default zero
% 'multitaper'        -- if true uses the multitaper method to estimate
%                        cross range spectrum
% OUTPUT
% 'P'                 -- voltage of output SAR image on meshgrid point of

% Here we have some hard coded parameters to handle some memory issues
cpisPerStep = 600;
gridPointsPerStep = 1E4;
C = 299792458;
% OSR is used to define how much oversampling is used in integration
% larger more accurate
% larger more slower
OSR = 16;
% defaults for various parameters the user can modify
cuda = false;
windowRangeType = @rectwin;
windowRangeParameters = [];
windowCrossRangeType = @rectwin;
windowCrossRangeParameters = [];
timewindow = [];
timeoverlap = 0;
multitaper  = false;
% get length of files
% read file header and get info on run
% Fs sample rate in Hz
% Fc carrier frequency in Hz
[samples_received, Fs, Fc] = testread(sprintf('%s',fileNameBase),-1);
% find min and max delay of saved data
pulseLength = length(pulseReplica)/Fs;
[systemDelay(2), systemDelay(1)] = finddelays( track,geometry,C,pulseLength);
%find samples in each pulse of processed data (based on how simulator
%works)
L = floor((systemDelay(2)-systemDelay(1))*Fs);
if(L<2*length(pulseReplica))
    L = 2*length(pulseReplica );
end
L = 2.^(nextpow2(L));

%block length is number of stored sample per pulse
blockLength =L;%-floor(length(pulseReplica)/2);
%calculate total number of pulses in file
totalBlocks = floor(samples_received/L)-1;
time = totalBlocks/prf;
for ii = 1:length(varargin)
    if(ischar(varargin{ii}))
        switch lower(varargin{ii})
            case {'cuda'}
                cuda = varargin{ii+1};
                if cuda ~=true&cuda~=false
                    
                    ME = MException('bpm:cuda', ...
                        sprintf('cuda must be set to true or false'));
                    throw(ME);
                end
            case {'windowrange'}
                wx= varargin{ii+1};
                if(iscell(wx))
                    windowRangeType = wx{1};
                    windowRangeParameters = wx(2:end);
                else
                    windowRangeType = wx;
                    windowRangeParameters = [];
                end
            case {'windowcrossrange'}
                wx= varargin{ii+1};
                if(iscell(wx))
                    windowCrossRangeType = wx{1};
                    windowCrossRangeParameters = wx(2:end);
                else
                    windowCrossRangeType = wx;
                    windowCrossRangeParameters = [];
                end
            case {'timewindow'}
                timewindow = varargin{ii+1};
                if timewindow <0||timewindow >time
                    
                    ME = MException('bpm:timewindow', ...
                        sprintf('Window times inconsistent'));
                    throw(ME);
                end
            case{'timeoverlap'}
                timeoverlap = varargin{ii+1};
                if timeoverlap <0||timeoverlap >time
                    
                    ME = MException('bpm:timewindow', ...
                        sprintf('Window times inconsistent'));
                    throw(ME);
                end
            case{'multitaper'}
                multitaper = varargin{ii+1};
                if multitaper  ~=true&multitaper ~=false
                    
                    ME = MException('bpm:multitaper', ...
                        sprintf('multitaper must be true or false'));
                    throw(ME);
                end
        end
    end
end

if(cuda)
    cast = @gpuArray;
    casttype = @single;
    inverse_cast = @gather;
    
else
    cast = @double;
    casttype = @double;
    inverse_cast = @deal;
end


[X,Y,Z] = meshgrid(xi,yi,zi);



if(~isempty(timewindow)&&~multitaper)
    %get lengths of overlapping windows and pulse to step between windows
    
    cpisInWindow = floor(timewindow * prf);
    cpisToAdvance = floor((timewindow-timeoverlap) * prf);
    
    windowBlocks = ceil(totalBlocks/cpisToAdvance);
    cpisToAdvance = totalBlocks/windowBlocks;
    fprintf(1,'***window size: %d advancing by %d total blocks: %d\n', cpisInWindow,cpisToAdvance,windowBlocks);
    if(cpisInWindow<3)
        
        ME = MException('bpm:timewindow', ...
            sprintf('cpisInWindow is too low'));
        throw(ME);
    end
    if(cpisToAdvance<3)
        ME = MException('bpm:timewindow', ...
            sprintf('cpisInWindow is too low'));
        throw(ME);
        
    end
    
    
else
    cpisInWindow = totalBlocks;
    windowBlocks = 1;
    cpisToAdvance = cpisInWindow;
end

%%construct fast time matched filter
%apply range sidelobe weighting and normalize
if(isempty(windowRangeParameters))
    weighting  = windowRangeType(length(pulseReplica));
else
    if(length(windowRangeParameters)==1)
        weighting  = windowRangeType(length(pulseReplica),windowRangeParameters{1});
    else
        weighting  = windowRangeType(length(pulseReplica),windowRangeParameters{1},windowRangeParameters{2});
    end
end

filterz  = (pulseReplica).* weighting;
filterz =  (filterz./norm(filterz));
%put in wrap around order and fft
h = zeros(1,L);
h(1:length(filterz)) = conj(filterz(end:-1:1));
h = circshift(h,[-length(filterz)/2,0]);
% h =zeros(size(h));
% h(1) = 1;
h = cast(casttype((fft(h.'))));

%make cross range window
if(~multitaper)
    if(isempty(windowCrossRangeParameters))
        weightingCross  = windowCrossRangeType(cpisInWindow);
    else
        if(length(windowCrossRangeParameters)==1)
            weightingCross  = windowCrossRangeType(cpisInWindow,windowCrossRangeParameters{1});
        else
            weightingCross  = windowCrossRangeType(cpisInWindow,windowCrossRangeParameters{1},windowCrossRangeParameters{2});
        end
    end
    weightingCross = weightingCross/norm(weightingCross);
else
    [E,V] = dpss(cpisInWindow,4);
    weightingCross = sqrt(E/sqrt(cpisInWindow));
end
if(~multitaper)
    P = casttype(zeros(size(X)));
    Ps = P;
else
    P = casttype(zeros(size(X,1),size(X,2),length(V)));
    
end
%phase term
K = cast(1i * 4 * pi  * Fc);
% we loop through in blocks to avoid memory/computation problems
% you can adjust cpisPerStep to work with the memory you have available.

% compute fast time range step size
dr = median(diff(ranges));
gridPointsCartesian = -cast(casttype([X(:)  Y(:)  Z(:)]));
%time between pulses
cpiTimeStep = round(Fs*1/prf)/Fs;

% total number of blocks we're dividing problem into (avoid memory problem)
megaBlocks = ceil(totalBlocks/cpisPerStep);
%how many cpis do we have left...

for window = 1:windowBlocks
    windowCpiStart = (window-1)*cpisToAdvance;
    pulsesRemaining = cpisInWindow;
    megaBlocks = ceil(pulsesRemaining /cpisPerStep);
    for mega = 1:megaBlocks
        cpiInMegaBlock = min([cpisPerStep  pulsesRemaining]);
        %read in block of cpisPerStep cpis
        offset = blockLength*(mega-1)*cpisPerStep+windowCpiStart;
        [dataReceived] = testread(sprintf('%s',fileNameBase),offset * 8,blockLength*cpiInMegaBlock);
        dataReceived = cast(reshape(dataReceived,blockLength,[]));
       % dataReceived = dataReceived + 1/sqrt(2)*(randn(size(dataReceived)) + 1i*randn(size(dataReceived)) );
        % perform fast time pulse compression
        p= ifft(bsxfun(@times,h,fft(dataReceived,L)));
       % p = p(length(filterz)/2+1:end,:);
      % p = p(length(filterz)/2:end,:);
        blockIndex = [0:(size(p,2)-1)]+(mega-1)*cpisPerStep;
        blockTime=  blockIndex*cpiTimeStep+(systemDelay(2)-systemDelay(1))/2 +systemDelay(1) ;
        %times in track
        t = track(4,:,1);
        %interpolate track of transmitter and receiver to times for pulses
        x_tr = cast(casttype(interp1(t ,track(1,:,1),blockTime,'linear','extrap')));
        y_tr = cast(casttype(interp1(t ,track(2,:,1),blockTime,'linear','extrap')));
        z_tr = cast(casttype(interp1(t ,track(3,:,1),blockTime,'linear','extrap')));
        
        x_r = cast(casttype(interp1(t ,track(1,:,2),blockTime,'linear','extrap')));
        y_r = cast(casttype(interp1(t ,track(2,:,2),blockTime,'linear','extrap')));
        z_r = cast(casttype(interp1(t ,track(3,:,2),blockTime,'linear','extrap')));
        
        % find range at given cpi for all points in grid
        % range is range from transmitter to grid point and then
        % grid point to receiver
        gridBlocks = ceil(size(gridPointsCartesian,1)/gridPointsPerStep);
        if(~multitaper)
            iss =  cast(casttype(zeros(1,size(gridPointsCartesian,1))));
        else
            iss =  cast(casttype(zeros(size(gridPointsCartesian,1),length(V))));
        end
        for gb = 1:gridBlocks
            indexGridBlock = [1:gridPointsPerStep]+(gb-1)*gridPointsPerStep;
            %trim off extra index points
            indexGridBlock(indexGridBlock>size(gridPointsCartesian,1)) = [];
            r=(sqrt( bsxfun(@plus,x_tr,gridPointsCartesian(indexGridBlock,1)).^2 +...
                bsxfun(@plus,y_tr,gridPointsCartesian(indexGridBlock,2)).^2 +...
                bsxfun(@plus,z_tr,gridPointsCartesian(indexGridBlock,3)).^2)+...
                sqrt( bsxfun(@plus,x_r, gridPointsCartesian(indexGridBlock,1)).^2 +...
                bsxfun(@plus,y_r, gridPointsCartesian(indexGridBlock,2)).^2 +...
                bsxfun(@plus,z_r, gridPointsCartesian(indexGridBlock,3)).^2))'*.5;
            ts = repmat(cast(casttype(1:size(r,1))),size(r,2),1)';
            %**trace integration path through data for given grid point
            % (normalized to integer grid)
            ri = (r-ranges(1))/dr;
            %*cut out region of interest for data based on evaluation grid points
            %find max and min data regions
%             [~,ind] = max(abs(p(:,1)));
%             ind = [-2:2]+ind;
%             pz = polyfit(ind,abs(p(ind,1))',2);
%             range = roots(polyder(pz));
            maxri = max(ri(:));
            minri = min(ri(:));
            %add ten sample buffer to each end
            region = [floor(minri)-100:ceil(maxri)+100];
            %cut region out of matrix
            pregion = double(gather(p(region,:)));
            %upsample to desired OSR prior to integration
            pregion = cast(casttype(resample(pregion ,OSR,1)));
            %make ranges of clipped data region
            rangeUp  = linspace(ranges(region(1)),ranges(region(end)),size(pregion,1));
            drx = (1-1/OSR)*dr/2;
            %make time index of data
            tx = 1:size(pregion,2);
            %range integration path
            %find values in image on given grid point paths
            s =  interp2(tx,rangeUp,pregion,ts,r-drx,'linear',0);
            %multiply by phase term
            s = s.*exp(K * (r/C));
            %apply window
            if(~multitaper)
                s = bsxfun(@times,s,weightingCross(blockIndex+1));
                %do intergration along grid point paths
                iss(indexGridBlock) = sum(s,1);
            else
                for mt = 1:size(weightingCross,2)
                  
                    %do intergration along grid point paths
                    iss(indexGridBlock,mt) = sum(bsxfun(@times,s,weightingCross(blockIndex+1,mt)),1);
                end
            end
            
        end
        
        %reshape and add to previous batch of pulses
        if(~multitaper)
            P =P+inverse_cast(reshape(iss,size(X,1),size(X,2),size(X,3)));
        else
            for mt = 1:size(weightingCross,2)
                P(:,:,mt) =P(:,:,mt)+inverse_cast(reshape(iss(:,mt),size(X,1),size(X,2),size(X,3)));
            end
        end
        pulsesRemaining = pulsesRemaining-cpisPerStep;
        fprintf(1,'Percent done: %3.2f\n',min([100 100 * mega/megaBlocks]));
    end
    %divide by total number of pulses
    P = gather(P);
    if(~multitaper&&windowBlocks>1)
        Ps = Ps + abs(P).^2;
        P = casttype(zeros(size(X)));
    end
end
if(multitaper)
    Ps = P;
    for mt = 1:2%size(weightingCross,2)
        P = abs(Ps(:,:,1)).^2;
    end
    P = P/mt;
    P = sqrt(P);
elseif(windowBlocks>1)
    P = sqrt(Ps/windowBlocks);
end
%end
% sr = real(s);
%         srgto = sr;
%         srlto = sr;
%         srgto(srgto<0) = 0;
%         srlto(srlto>0) = 0;
%         si = imag(s);
%         sigto = si;
%         silto = si;
%         sigto(sigto<0) = 0;
%         silto(silto>0) = 0;
%         iss(indexGridBlock) = sum(srgto,1)+sum(srlto,1)+1i*...
%                               (sum(silto,1)+sum(sigto,1));
