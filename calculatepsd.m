
function [spectrum, freqs, time] = calculatepsd(data,Fs,N,overlap,window,real)
start_index_inc =  floor(N * (1-overlap));
K = (floor((size(data,1)/(start_index_inc))));
if((K*start_index_inc+N) >size(data,1))
    K = (floor(((size(data,1)-2*N)/(start_index_inc))));
end
if (K<1) 
    errormessage = sprintf('Insufficent samples , %d, to process run given fft length of %d.  Increase cut length or decrease fft size',...
        size(data,1),N);
    ME = MException('VerifyOutput:OutOfBounds', ...
        errormessage);
    throw(ME);
end
if (length(window)~=N) 
    errormessage = sprintf('Incorrect samples , %d, in window.  Window size must equal FFT size, %d',...
        length(window),N);
    ME = MException('VerifyOutput:OutOfBounds', ...
        errormessage);
    throw(ME);
end
start_index = 1;
% create data matrix
window = window * norm(ones(N,1))/norm(window);
window = repmat(window,1,size(data,2));
if(real)
    xf = single(zeros(N/2,size(data,2),K));
    for ii = 1:K
        f =  fft(window.*data(start_index:N+start_index-1,:))/N;
        xf(:,:,ii) = f(2:end/2+1,:)'.';
        time(ii) = ((ii-1)) * start_index_inc/Fs +N/2/Fs;
        start_index = start_index + start_index_inc;

    end
    freqs = Fs/2*linspace(0,1,NFFT/2);
else
    xf = single(zeros(N,size(data,2),K));

    for ii = 1:K-1
        
        f =  fftshift(fft(window.*data(start_index:(N+start_index-1),:))/N,1);
        xf(:,:,ii) = f'.';
        time(ii) = ((ii-1)) * start_index_inc/Fs+N/2/Fs;
        start_index = start_index + start_index_inc;
    end
    freqs  = Fs * [-N/2:(N/2-1)]/N;

    
end
spectrum = xf;