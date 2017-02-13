function [w power] = element_stap_2(rsp,steering_vector,Npulses,snapshots,guard_bandwidth,points_kept)
Wd = fftshift(dftmtx(Npulses),2);
rsp = permute(rsp,[1 3 2]);
L = ceil(snapshots/2);

count = 1;
Nsensors = size(rsp,2);
for range_gate =100:size(rsp,1)
    start_left = max([range_gate-(L+guard_bandwidth) 1 ]);
    end_left = max([range_gate - guard_bandwidth-1 1]);
    start_right = min([range_gate + guard_bandwidth+1 size(rsp,1)-1]);
    end_right = min([range_gate+(L+guard_bandwidth) size(rsp,1)-1]);
    s = (rsp([start_left:end_left],:,:));
    s2 = abs(s(:,:,floor(size(s,3)/2)));
    [mv iii]=  sort(max(s2,[],2));
    s = s(iii(end:-1:max([end-points_kept 1])),:,:);
    % s = s(iii(end),:,:);
    s = permute(reshape(s(:),size(s,1),Nsensors,Npulses,[]),[1 4 2 3]);
    s = reshape(s(:),[],Nsensors *Npulses);
    sl = s;
    s = (rsp([start_right:end_right],:,:));
    s2 = abs(s(:,:,floor(size(s,3)/2)));
    [mv iii]=  sort(max(s2,[],2));
    s = s(iii(end:-1:max([end-points_kept 1])),:,:);
    %s = s(iii(end),:,:);
    s = permute(reshape(s(:),size(s,1),Nsensors,Npulses,[]),[1 4 2 3]);
    s = reshape(s(:),[],Nsensors * Npulses);
    sr = s;
    s = [sl;sr];
    R = (s' * s)/size(s,1);
    R = R+1 * eye(size(R));
   % R = R + eye(size(R));
% 
    Ri = inv(R);
   
    s = (rsp(range_gate,:,:));
    s = permute(reshape(s(:),size(s,1),Nsensors,Npulses,[]),[1 4 2 3]);
    s = reshape(s(:),[],Nsensors * Npulses).';
    for ii = 1:size(Wd,2)
        steering_vector_doppler = Wd(:,ii);
        sv = kron(steering_vector_doppler,steering_vector);%
        w = Ri * sv;
        
        power(count,ii) = w' * s;
    end
    count = count + 1
end
