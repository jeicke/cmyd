% pulse_compression filtering of radar data by replica
% example
% [D] = pulse_compression(rsp,radar,'cheb',70);
function [D] = pulse_compression(rsp,pulse_replica,compression_weighting,beta)
D = double(zeros(size(rsp)));
switch lower(compression_weighting)
    case('none')
        weighting  = rectwin(length(pulse_replica));
    case('hanning')
        weighting = (0.54+0.46*cos(2*pi*[(-length(pulse_replica)/2)+0.5:(length(pulse_replica)/2)-0.5]/length(pulse_replica)))';
    case('cheb')
        weighting  =  chebwin(length(pulse_replica),beta) ;
       %weighting = chebwgt(length(pulse_replica),beta);
    otherwise
        weighting = (0.54+0.46*cos(2*pi*[(-length(pulse_replica)/2)+0.5:(length(pulse_replica)/2)-0.5]/length(pulse_replica)))';
end
weighting =double(weighting);
pulse_replica =  double(pulse_replica./norm(pulse_replica));
if(length(pulse_replica)>1)
    filterz  = conj(pulse_replica).* weighting;
    for l=1:size(rsp,3)

        %(:,:,l)  = fftfilt(filterz,rsp(:,:,l));
        D(:,:,l)  = filter(filterz,1,rsp(:,:,l));
    end
else
    D = rsp;
end
