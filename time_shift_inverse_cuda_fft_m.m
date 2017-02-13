% N is fft size
% v is speed of light
function [D P K] = time_shift_inverse_cuda_fft_m(Fs,N,V,pulse_replica,ranges,operating_frequency,P,K)



if(isempty(K))
    f = ifftshift((-N/2:N/2-1)'/N*Fs);
    K = gpuArray(single(-1i * 4 * pi  * (f - operating_frequency)));
end
%H = H.*repmat((exp(-1i * operating_frequency * 4 * pi * ranges/V)),length(f),[]);
if(isempty(P))
    P = fft(pulse_replica ,N);
    P = gpuArray(single(repmat(P,size(ranges,2),[]).'));
end
D = P .* exp(K *gpuArray(single(ranges/V )));

%D = D2(1:N2,:);
