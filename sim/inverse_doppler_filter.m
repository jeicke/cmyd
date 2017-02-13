% doppler filtering of radar data by replica
% example
% [D] = inverse_doppler_filter(rsp,128);
function [D] = inverse_doppler_filter(rsp)


D= single(zeros(size(rsp)));

for k=1:size(rsp,3)
    D(:,:,k)=ifft(squeeze(rsp(:,:,k).')).';
end
