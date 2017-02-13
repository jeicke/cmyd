function [Dbf, w] = beamform_stap(D,Vs,Vd,snapshots,guard_bandwidth,sigma,points_kept)
[w] = pre_doppler_stap_2(D,1:size(D,3),Vs,Vd,snapshots,guard_bandwidth,sigma,false,points_kept);
M = size(D,2)-size(w,1)/size(D,3)+1;
Dbf = single(zeros(size(D,1),M,size(w,2)));

for jj = 1:size(w,2)
    [D2] = apply_doppler_stap_filter_3(D,squeeze(w(:,jj,:)));
    %[D2 radar] = doppler_filter(D2,M,'cheb',70);
    Dbf(:,:,jj) = D2;
end
w2 = reshape(w(:),size(D,3),size(w,1)/size(D,3),[]);
w = single(w2);%single(w2(:,end:-1:1,:));


