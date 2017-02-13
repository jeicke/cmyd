% function beamforms data
function [bf] = beamform(rsp,channels,steering_vector)
rsp = rsp(:,:,channels);
D2 = permute(rsp,[3 1 2]);
bf = single(zeros(size(rsp,2),size(steering_vector,2),size(rsp,1)));
for ii = 1:size(rsp,2)
    bf(ii,:,:) = steering_vector.' * D2(:,:,ii);
end

bf = permute(bf,[3 1 2]);