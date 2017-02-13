function [P,sinr,doppler_bins] =stap( fileNameBase,prf,pulsesincpi,pulseReplica,ranges,track,geometry,Vs,Vd,snapshots,guardband,velocity,Fc,varargin)
% BPM Performs SAR time domain correlation
% Stands for:
% Back Propogation Method (BPM)
% perfoms backpropogation of SAR data
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
% 'timeoverlap'       -- length of time in seconds that segments overlap
%                        default zero
% 'multitaper'        -- if true uses the multitaper method to estimate
%                        cross range spectrum
% OUTPUT
% 'P'                 -- voltage of output SAR image on meshgrid point of
% the SPEED of light
C = 299792458;
P=[];
% defaults for various parameters the user can modify
cuda = false;
windowRangeType = @rectwin;
windowCrossRangeType = @rectwin;
windowRangeParameters = [];
windowCrossRangeParameters = [];
method = 'STAPFULL';
dopplerwarping = [];
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
            case {'windowdoppler'}
                wx= varargin{ii+1};
                if(iscell(wx))
                    windowCrossRangeType = wx{1};
                    windowCrossRangeParameters = wx(2:end);
                else
                    windowCrossRangeType = wx;
                    windowCrossRangeParameters = [];
                end
             case {'dopplerwarping'}
                dopplerwarping= varargin{ii+1};
                
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
    weightingCross  = windowCrossRangeType(cpisToAdvance);
else
    if(length(windowCrossRangeParameters)==1)
        weightingCross  = windowCrossRangeType(cpisToAdvance,windowCrossRangeParameters{1});
    else
        weightingCross  = windowCrossRangeType(cpisToAdvance,windowCrossRangeParameters{1},windowCrossRangeParameters{2});
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
%make composity steering vector
BOM = kron(Vs,Vd);
opts.POSDEF = true; opts.SYM = true;
for mega = 1%:megaBlocks
    cpiInMegaBlock = min([cpisToAdvance  pulsesRemaining]);
    %read in block of cpisPerStep cpis
    offset = blockLength*(mega-1)*cpisToAdvance;
    [dataReceived] = testread(sprintf('%s',fileNameBase),offset * 8,blockLength*cpiInMegaBlock,blockLength);
    for channel = 1:size(dataReceived,2)
        datas= cast(squeeze(dataReceived(:,channel,:)));%cast(reshape(dataReceived(:,channel),blockLength,[]));
        % dataReceived = dataReceived + 1/sqrt(2)*(randn(size(dataReceived)) + 1i*randn(size(dataReceived)) );
        % perform fast time pulse compression
        p(:,:,channel)= ifft(bsxfun(@times,h,fft(datas,L)));
        
        % p = p(length(filterz)/2:end,:);
        
        %times in track
        
        
    end
    dim = size(p);
    %matrix multiply to do range-Doppler process
    switch method
        case 'CONV'
            p = reshape(p,dim(1),dim(2)*dim(3))*BOM;
            p = reshape(p,dim(1),dim(2),size(Vs,2));
        case 'STAPFULL'
            
            for ii = 295
               
                % estimate spatial-doppler correlation matrix for each range using N snapshots
                % around given range value with a given guard band
                startLeft = max([ii-(snapshots+guardband) 1 ]);
                endLeft = max([ii - guardband-1 1]);
                startRight = min([ii + guardband+1 size(p,1)-1]);
                endRight = min([ii+(snapshots+guardband) size(p,1)-1]);
                
                % create interference plus noise covariance matrix using size a 3
                % weight doppler filter
                
                s = conj(p([startLeft:endLeft startRight:endRight],:,:));
                dims = size(s);
                if(~isempty(dopplerwarping))
                    for rg= 1:numel(dopplerwarping)
                        p1 = track(1:3,1,1);
                        p2 = track(1:3,1,2);
                        tp = [dopplerwarping(rg) 0 0];
                        r(rg) = norm(p1-p2)+norm(p2+tp);
                        a =1 ;
                    end
                     [clutter_doppler, EL] = create_doppler_bi(90*pi/180,velocity',dopplerwarping,Fc,track(3,1,2)/1000);
                      dv = interp1(r,clutter_doppler,ranges([startLeft:endLeft startRight:endRight]),'linear');
                      indx = ~isnan(dv)
                       s = s(indx,:,:);
                      %[clutter_dopplertg, EL] = create_doppler_bi(0,velocity',ranges(ii),Fc,track(3,1,2)/1000);
                     % correction = clutter_dopplertg-clutter_doppler;
                     pulse=((0:(size(p,2)-1))')/prf;
                     dv = dv(indx);
                     target = mean(dv);
                     cValue = (target -dv);
                     
                     for dvs = 1:size(s,1)
                         s2(dvs,:,:) =bsxfun(@times,squeeze(s(dvs,:,:)),exp(-1i*2*pi*pulse*cValue(dvs)));
                     end
                     a=1;
                end
                dims = size(s);
                s = reshape(s(:),dims(1),dims(2)*dims(3));
                
                R = (s' * s)/size(s,1);
            end
            % solve for weights
           
            R = (s' * s)/size(s,1);
            
            pa = [];
            doppler_bins = [-prf/2:2:prf/2];
            %[sinr,w] = sinrR(R,pulsesincpi,7,doppler_bins ,Vs(:,1),prf);
            [sinr{1}] = sinrR(R,pulsesincpi,size(Vs,1),doppler_bins ,Vs(:,1),prf);
            [sinr{2}] = sinrRc(R,pulsesincpi,size(Vs,1),doppler_bins ,Vs(:,1),prf);
            [sinr{3}] = sinrRpdm2(R,pulsesincpi,size(Vs,1),doppler_bins ,Vs(:,1),prf,1);
            [sinr{4}] = sinrRpdm2(R,pulsesincpi,size(Vs,1),doppler_bins ,Vs(:,1),prf,2);
           % pa = w'*s';
            a = 1;
        otherwise
    end
    p = pa;
    
end


% pf = fft(bsxfun(@times,p,weightingCross(1:cpiInMegaBlock).').',cpiInMegaBlock);

P(:,:,:,mega) = p;
blockIndex = [0:(size(p,2)-1)]+(mega-1)*cpisToAdvance;
blockTime=  blockIndex*cpiTimeStep+(systemDelay(2)-systemDelay(1))/2 +systemDelay(1) ;
pulsesRemaining = pulsesRemaining-cpisToAdvance;
fprintf(1,'Percent done: %3.2f\n',min([100 100 * mega/megaBlocks]));
end
% or ta = 1:size(s,1)
%                          test = squeeze(s(ta,:,:));
%                          x = test*Vs(:,1);
%                          z(:,ta) = fftshift(20*log10(abs((ifft(x)))));
%                          
%                          [h,i(ta)] = max(z(:,ta));
%                      end
%                      ic =33;% mode(i);
%                      a = 1;
%                      correction = (i-ic);
% 
%                      for ta = 1:size(s,1)
%                          %dv = exp(1i*2*pi*pulse*correction(ta) );
%                          test = squeeze(s(ta,:,:));
%                          for ii = 1:size(test,2)
%                              test(:,ii) = fftshift(ifft(test(:,ii)));
%                              test(:,ii) = circshift(test(:,ii),-correction(ta));
%                              test(:,ii) = fft(ifftshift(test(:,ii)));
%                          end
%                          %correct = bsxfun(@times,dv,test);
%                          s(ta,:,:) = test;
%                          z1(:,ta) = fftshift(20*log10(abs((ifft( test*Vs(:,1))))));
%                      end
%                      a = 1;