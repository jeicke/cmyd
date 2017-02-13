function [D] = apply_doppler_stap_filter_3(rsp,w)




M = size(rsp,2)-size(w,1)/size(rsp,3)+1;

w2 = reshape(w(:),size(rsp,3),size(w,1)/size(rsp,3),[]);
w2 = single(w2(:,end:-1:1,:));
%w is sensor x doppler window vector in order Sensor1(dw1 dw2 dw3...)
%Sensor2(dw1 dw2 dw3) etc
d = single(zeros(size(rsp)));

for channel = 1:size(rsp,3)
     data = squeeze(rsp(:,:,channel));
     w3 = squeeze(w2(channel,:,:));
    %d(:,:,channel) = fftfilt(w3(:,:),data(:,:).').';
    for cpi = 1:size(rsp,1)
        
       d(cpi,:,channel) = filter(w3(:,cpi),1,data(cpi,:)).';
       
%         d2(cpi,:,channel) = 0;
%         for p = 1:M
%             d2(cpi,p+2,channel) = d2(cpi,p+2,channel)+w3(:,cpi).'*data(cpi,(1:3)+(p-1)).';
%             %filter(w3(:,cpi),1,data(cpi,:)).';
%         end
    end
  

end

 D = sum(d(:,3:end,:),3);
