function [P, blockTimes] =processradarblock( fileNameBase,pulse_replica,systemDelay)
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

[samples_received, Fs, Fc, channels] = testread(sprintf('%s',fileNameBase),-1);

L = floor((systemDelay(2)-systemDelay(1))*Fs)+floor(length(pulse_replica)/2);
if(L<2*length(pulse_replica))
    L = 2*length(pulse_replica );
end
L = 2.^(nextpow2(L));


blockLength =L-floor(length(pulse_replica)/2);
totalBlocks = floor(samples_received/L)-1;
% 
blocks = 0;
state = [];



totalBlocks= floor(totalBlocks);
p = zeros(L,1,floor(totalBlocks));
blockTimes = zeros(1,totalBlocks);
%weighting  =  hannwindow(length(pulse_replica)) .';
weighting  = chebwin(length(pulse_replica),70); 
pulse_replica =  double(pulse_replica./norm(pulse_replica));
filterz  = (pulse_replica).* weighting;
h = zeros(1,L);
h(1:length(filterz)) = filterz;
h = circshift(h.',-length(filterz)/2);
h = conj(fft(h));
while blocks<totalBlocks
   
    offset = blockLength*blocks+1;

    tStart = offset/Fs;
    tEnd = (offset + blockLength-1)/Fs;
    [data_received] = testread(sprintf('%s',fileNameBase),...
        offset * 8,blockLength );

  %  data_received = data_received  + 1/sqrt(2) * (randn(size(data_received ))+ 1i * randn(size(data_received )));
    times = (tStart:1/Fs:tEnd);
    state.create = true;

    D = data_received;
    
    
    D = ifft(fft(D,L).'.*h);
    %toc;
    %D = [zeros(1,length(pulse_replica)) data_received zeros(1,length(pulse_replica))];
    %[D ] = pulsecompression(D,pulse_replica,'cheb',50).';
    
    %D = D(length(pulse_replica):end-length(pulse_replica));
    %D = D(rangeGates,:);
    
    %D = D(length(pulse_replica):end-length(pulse_replica));
    D = D(:,:);
    
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

