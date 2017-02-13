% doppler filtering of radar data by replica
% example
% [D] = doppler_filter(radar,rsp,128,'cheb',70);
function [D  shade_doppler] = doppler_filter(rsp,fftsize,doppler_weighting,beta,shift_fft)
if(~exist('shift_fft','var'))
    shift_fft = true;
end
switch lower(doppler_weighting)
    case('none')
        shade_doppler  = rectwin(size(rsp,2));
    case('hanning')
        shade_doppler = (0.54+0.46*cos(2*pi*[(-size(rsp,2)/2)+0.5:(size(rsp,2)/2)-0.5]/size(rsp,2)))';
    case('cheb')
        shade_doppler =  chebwin(size(rsp,2),beta) ;
        %shade_doppler = chebwgt(size(rsp,2),beta);
    otherwise
        shade_doppler = (0.54+0.46*cos(2*pi*[(-size(rsp,2)/2)+0.5:(size(rsp,2)/2)-0.5]/size(rsp,2)))';
end
shade_doppler = single(shade_doppler);

D= single(zeros(size(rsp,1),fftsize,size(rsp,3)));
if(shift_fft )
    for k=1:size(rsp,3)
        D(:,:,k)=fftshift(fft(squeeze(rsp(:,:,k).').*repmat(shade_doppler,1,size(rsp,1)) ,fftsize),1).';
    end
else
    for k=1:size(rsp,3)
        D(:,:,k)=(fft(squeeze(rsp(:,:,k).').*repmat(shade_doppler,1,size(rsp,1)) ,fftsize)).';
    end
end
