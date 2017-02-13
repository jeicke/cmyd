function [P blockTimes] =processradar( fileNameBase,prf,cpi,pulse_replica,rangeGates,time,fftsize)
%processes radar data from simulation
% fileNameBase        -- base filename for data to process
% 'prf'               -- prf of radar pulses
% 'nPulses'           -- number of pulses to use in cpi
% 'pulseBandwidth'    --
% 'pulseLength'       -- pulse length in seconds
% 'time'              -- [start end] time to process
% 'cpi'               -- Number of pulses in cpi
cuda = false;
if(cuda)
    cast = @GPUsingle;
    inverse_cast = @double;
else
    cast = @double;
    inverse_cast = @double;
end
parallel = false;


P = [];
% get length of files
if(~exist('C','var')||isempty(C))
    C = 299792458;
end
if(~exist('coeff','var')||isempty(coeff))
    coeff = false;
end
%track = load(sprintf('%s_%s',fileNameBase,'track.mat'));

[samples_received Fs Fc channels] = testread(sprintf('%s',fileNameBase),-1);

% make replica pulse

plength = round(Fs/ prf);


% track = load(sprintf('%s_%s',fileNameBase,'track.mat'));
% 
if(exist('time','var')&&~isempty(time))
    start_sample = floor(time(1) * Fs);
    if(start_sample>samples_received)
        fprintf(2,'ERROR: TIME IS GREATER THAN TIME IN FILE');
        return;
        
    end
    end_sample = ceil(time(2) * Fs);
    if(end_sample>samples_received)
        end_sample = samples_received;
    end
else
    start_sample = 0;
    end_sample = samples_received;
end
% 
 blockLength =cpi * plength;
 
 totalBlocks = (end_sample-start_sample)/blockLength;
% 
blocks = 0;
state = [];
state.maxSystemDelay = blockLength  *4/Fs;
% 
if(~exist('rangeGates','var')||isempty(rangeGates))
    rangeGates = 1:plength;
end
if( parallel)
    matlabpool open
end
data_received(1) = 1;
first = true;
max_system_delay  = floor(state.maxSystemDelay * Fs)/Fs;
if(~exist('fftsize','var')|isempty(fftsize))
    fftsize = cpi;
end
totalBlocks = floor(totalBlocks);
p = zeros(length(rangeGates),fftsize ,totalBlocks);
blockTimes = zeros(1,totalBlocks);
weighting  =  hannwindow(length(pulse_replica)) .';
pulse_replica =  double(pulse_replica./norm(pulse_replica));
filterz  = conj(pulse_replica).* weighting;

while blocks<totalBlocks
    
    tStart = (blocks * blockLength + start_sample)/Fs;
    
    tEnd = (start_sample + blocks * blockLength + blockLength - 1)/Fs;
    blockTimes(blocks+1) = tStart;
    
    if(tEnd * Fs>end_sample)
        tEnd =end_sample/Fs;
        blockLength = ceil((tEnd-tStart) * Fs);
    end


    offset = round(tStart * Fs);

    tStart = offset/Fs;
    tEnd = (offset + blockLength-1)/Fs;
    [data_received] = testread(sprintf('%s',fileNameBase),...
        offset * 8,blockLength );

    %data_received = data_received  + 1/sqrt(2) * (randn(size(data_received ))+ 1i * randn(size(data_received )));
    times = (tStart:1/Fs:tEnd);
    state.create = true;
    
    %datapad = [ data_received  zeros(1,ceil((length( data_received )/ plength)/cpi) *  plength*cpi-length( data_received ))];
    %D = reshape(datapad,plength,cpi,[]);
   
   % tic;
    D = [zeros(1,length(pulse_replica)) data_received];

    if(first)
         n = nextpow2(length(D));
        h = fft(filterz,2^n);
    end
    D = ifft(fft(D,2^n).'.*h);
    %toc;
    %D = [zeros(1,length(pulse_replica)) data_received zeros(1,length(pulse_replica))];
    %[D ] = pulsecompression(D,pulse_replica,'cheb',50).';
    
    %D = D(length(pulse_replica):end-length(pulse_replica));
    %D = D(rangeGates,:);
    
    %D = D(length(pulse_replica):end-length(pulse_replica));
    D = D(rangeGates,:);
    
    blocks = blocks + 1;
    p(:,:,blocks) = D(:,:);
    fprintf(1,'Percent done: %3.2f\n',min([100 100 * blocks/totalBlocks]));
   % p = 20*log10(abs(P));
        first = false;
end
if( parallel)
    matlabpool close
end
P = inverse_cast(p);

% save(sprintf('%s_%s',fileNameBase,'correlation.mat'),'P');
% end

