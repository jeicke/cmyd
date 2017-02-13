function [H] = delayFilterFreq(K,ranges,sigma,Kfc,cast,Tz)
step_size =1000;
%get index
start_index = 1;
end_index = step_size;
if(end_index>size(ranges,2))
    end_index = size(ranges,2);
end
index = start_index:end_index;
% if(~(size(ranges,1)==1))
%     K = repmat(K,1,size(ranges,1));
% end

H = cast(single(zeros(size(K,1),size(ranges,1))));
while (start_index<=size(ranges,2))
    r = ranges(:,index);
    s = (sigma(:,index));
    M= (s.*exp(Kfc*r)).';
    r = r-Tz;
    if(size(ranges,1)==1)
        H = exp( K* r)*M;
    else
        for ii = 1:size(ranges,1)
            H(:,ii) = H (:,ii)+ exp( K* r(ii,:))*M(:,ii);
        end
    end
    % update index
    start_index = start_index + step_size;
    end_index = end_index + step_size;
    index = start_index:end_index;
    if(end_index>size(ranges,2))
        index = start_index:length(ranges);
    end
end