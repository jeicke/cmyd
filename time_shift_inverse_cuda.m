% N is fft size
% v is speed of light
function [D P] = time_shift_inverse_cuda(Fs,N,N2,V,pulse_replica,ranges,operating_frequency,P)



ranges = GPUsingle(ranges/V);

f = GPUsingle(ifftshift((-N/2:N/2-1)'/N*Fs));
H = exp(GPUsingle(-1i * 4 * pi  * (f + operating_frequency)) * ranges);


%H = H.*repmat((exp(-1i * operating_frequency * 4 * pi * ranges/V)),length(f),[]);
if(isempty(P))
    P = fft(pulse_replica ,N);
    P = repmat(P,size(H,2),[]).';
    P = GPUsingle(P);
end
D2 = ifft(P .* H);
D = D2(1:N2,:);


