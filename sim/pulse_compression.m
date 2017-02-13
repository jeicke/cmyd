% pulse_compression filtering of radar data by replica
% example
% [D] = pulse_compression(rsp,radar,'cheb',70);
function [D] = pulse_compression(rsp,pulse_replica,compression_weighting,beta)
D = single(zeros(size(rsp)));
switch lower(compression_weighting)
    case('none')
        weighting  = rectwin(length(pulse_replica));
    case('hanning')
        weighting = (0.54+0.46*cos(2*pi*[(-length(pulse_replica)/2)+0.5:(length(pulse_replica)/2)-0.5]/length(pulse_replica)))';
    case('cheb')
         %shade_doppler =  chebwin(size(rsp,2),beta) ;
        weighting = chebwgt(length(pulse_replica),beta);
    otherwise
        weighting = (0.54+0.46*cos(2*pi*[(-length(pulse_replica)/2)+0.5:(length(pulse_replica)/2)-0.5]/length(pulse_replica)))';
end
weighting =single(weighting);
pulse_replica =  single(pulse_replica./norm(pulse_replica));
if(length(pulse_replica)>1)
    filterz  = conj(pulse_replica).* weighting;
    for l=1:size(rsp,3)
        %for k = 1:size(rsp,2)
%              z = rsp(:,:,l);
%             z = [zeros(length(filterz),100); z; zeros(length(filterz),100)];
%             z2 = fftfilt(filterz,z);
%             D(:,:,l) = z2(length(filterz)+1:(end-length(filterz)),:);
            D(:,:,l)  = fftfilt(filterz,rsp(:,:,l));
            %D(:,:,l)  = filter(filterz,1,rsp(:,:,l));
        %end
    end
else
    D = rsp;
end
%for k = 1:size(rsp,2)
          %  z = rsp(:,:,l);
           % z = [zeros(length(filterz),100); z; zeros(length(filterz),100)];
            %D(:,:,l) = fftfilt(filterz,z);
            %D(:,:,l) = z2(length(filterz)+1:(end-length(filterz)),:);
            %D(:,:,l)  = filter(filterz,1,rsp(:,:,l));
        %end