function [P] =rps( fileNameBase,prf,pulsesincpi,pulseReplica,ranges,track,geometry,varargin)
% RPS Radar Processing System
% Does conventional radar processing
% INPUT
% 'fileNameBase'      -- base filename for data to process=
% 'prf'               -- Pulse repetition frequency in Hz
% 'pulseReplica'      -- IQ samples of transmitted fast time pulse
% 'geometry'          -- Geometry of receive array #elementsx3 (x,y,z)
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
% 'postdoppler'      -- Number of post doppler adaptive elements
% 'nofft'            --don't do fft
% OUTPUT
% 'P'                 -- voltage of output SAR image on meshgrid point of
% the SPEED of light
C = 299792458;
% defaults for various parameters the user can modify
cuda = false;
windowRangeType = @rectwin;
windowCrossRangeType = @rectwin;
windowRangeParameters = [];
windowCrossRangeParameters = [];
nofft = false;
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
pulses = 1;
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
            case {'windowdoppler'}
                wx= varargin{ii+1};
                if(iscell(wx))
                    windowCrossRangeType = wx{1};
                    windowCrossRangeParameters = wx(2:end);
                else
                    windowCrossRangeType = wx;
                    windowCrossRangeParameters = [];
                end
            case {'postdoppler'}
                pulses=varargin{ii+1};
            case {'nofft'}
                nofft=varargin{ii+1};
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


cpisInWindow = pulsesincpi;
windowBlocks = 1;
cpisToAdvance = pulsesincpi;

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
if(isempty(windowCrossRangeParameters))
    weightingCross  = windowCrossRangeType(cpisToAdvance-pulses+1);
else
    if(length(windowCrossRangeParameters)==1)
       weightingCross  = windowCrossRangeType(cpisToAdvance-pulses+1,windowCrossRangeParameters{1});
    else
        weightingCross  = windowCrossRangeType(cpisToAdvance-pulses+1,windowCrossRangeParameters{1},windowCrossRangeParameters{2});
    end
end
filterz  = (pulseReplica).* weighting;
filterz =  (filterz./norm(filterz));
%put in wrap around order and fft
h = zeros(1,L);
h(1:length(filterz)) = conj(filterz(end:-1:1));
%h = circshift(h,[-length(filterz)/2,0]);
% h =zeros(size(h));
% h(1) = 1;
h = cast(casttype((fft(h.'))));
%phase term
K = cast(1i * 4 * pi  * Fc);

%time between pulses
cpiTimeStep = round(Fs*1/prf)/Fs;

% total number of blocks we're dividing problem into (avoid memory problem)
megaBlocks = floor(totalBlocks/cpisToAdvance);
%how many cpis do we have left...
pulsesRemaining = totalBlocks;
for mega = 1:megaBlocks
    cpiInMegaBlock = min([cpisToAdvance  pulsesRemaining]);
    %read in block of cpisPerStep cpis
    offset = blockLength*(mega-1)*cpisToAdvance;
    [dataReceived] = testread(sprintf('%s',fileNameBase),offset * 8,blockLength*cpiInMegaBlock,blockLength);
    for channel = 1:size(dataReceived,2)
        datas= cast(squeeze(dataReceived(:,channel,:)));%cast(reshape(dataReceived(:,channel),blockLength,[]));
        % dataReceived = dataReceived + 1/sqrt(2)*(randn(size(dataReceived)) + 1i*randn(size(dataReceived)) );
        % perform fast time pulse compression
        p= ifft(bsxfun(@times,h,fft(datas,L)));
        if(pulses>1)
            for pulse = 1:max([1 pulses])
                endIndex = size(p,2)-(pulses-pulse);
                pf = fft(bsxfun(@times,p(:,pulse:endIndex),weightingCross.').',cpiInMegaBlock-pulses+1);
                P{pulse}(:,:,channel,mega) = pf;
            end
        else
            pulse = 1;
            if(~nofft)
                pf = fft(bsxfun(@times,p,weightingCross.').',cpiInMegaBlock);
            else
                pf = p;
            end
            P{pulse}(:,:,channel,mega) = pf;
        end
        
        % p = p(length(filterz)/2:end,:);
        blockIndex = [0:(size(p,2)-1)]+(mega-1)*cpisToAdvance;
        blockTime=  blockIndex*cpiTimeStep+(systemDelay(2)-systemDelay(1))/2 +systemDelay(1) ;
        %times in track
    end
    pulsesRemaining = pulsesRemaining-cpisToAdvance;
    fprintf(1,'Percent done: %3.2f\n',min([100 100 * mega/megaBlocks]));
end
