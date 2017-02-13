function [result] = pulsecompress(matchedfilter, data,L,P,cudaParams)
%DETECTTRANSMISSION Detects Radiomap known transmitted signal
% Input
% matchedfilter      -- transmitted signal in time domain
% data               -- received signal
% L                  -- Length of fft in filter
% P                  -- Number of filter taps
% Ouput
% result             -- Output pulse compressed data

%if first call pad with zeros

if(~exist('cudaParams','var')||isempty(cudaParams))
    cudaParams.arraytype =  @deal;
    cudaParams.datatype = @deal;
    cudaParams.arraytypeString = '';
    cudaParams.cuda = false;
end

% break input into fft sections using overlap-save
if(cudaParams.cuda&&~isa(data,'gpuArray'))
    data = cudaParams.arraytype(cudaParams.datatype(data(:)));
end

D = L-(P-1);
% pad data with zeros
NN = gpuArray(length(data));
if(cudaParams.cuda)
   data = [zeros(P-1,1,'single','gpuArray'); data(:)];
    PL = ceil(length(data)/D)*D+D;
    data = [data.' zeros(1,PL-length(data),'single','gpuArray')].';
    
else
    data = [zeros(1,P-1) data(:).'].';
    PL = ceil(length(data)/D)*D+D;
    data = [data.' zeros(1,PL-length(data))].';
end

endl = fix(NN/D-1)+1;
if(endl*D+L>length(data))
    endl = endl-1;
end
if(cudaParams.cuda)
    Xrr = zeros(L,endl+1,'single','gpuArray');
else
    Xrr = zeros(L,endl+1);
end
totalLength = 0;

% for r = 0:endl
%     
%     index = r*D + 1:r*D+L;
%     
%     xr = data(index);
%     Xrr(:,r+1) = xr;
%     totalLength = totalLength + D;
% end

index = bsxfun(@plus,gpuArray([0:L-1]'),gpuArray(1:D:(endl+1)*D));
Xrr = fft(data(index));
H =fft(matchedfilter,L);

%yj = ifft( Xrr.*repmat(H,1,size(Xrr,2)));
yj = (ifft( bsxfun(@times,Xrr,H)));
%index = fix(P/2):L-fix(P/2);
index = P:L;
y  = yj(index,:);
y = y(:);
result = y(1:NN);