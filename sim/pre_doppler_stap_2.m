function [w] = pre_doppler_stap2(rsp,channels,steering_vector,steering_vector_doppler,snapshots,guard_bandwidth,sigma,max_power,points_kept)
count = 1;
rsp = rsp(:,:,channels);
% create kronicker product steering vector
C = zeros(size(rsp,3) * 3,size(steering_vector_doppler,2) * size(steering_vector,2));
L = ceil(snapshots/2);
for jj = 1:size(steering_vector_doppler,2)
    for ii = 1:size(steering_vector,2)



        C(:,count) = colkron(steering_vector_doppler(:,jj),steering_vector(:,ii));
        count = count + 1;
    end
end
%rsp = reshape(rsp(:),size(steering_vector_doppler,1),[]);
if(~exist('max_power','var'))
    max_power = false;
end
opts.POSDEF = true; opts.SYM = true;
snapshots = floor(size(rsp,2)/size(steering_vector_doppler,1)) * size(steering_vector_doppler,1); 
w = zeros(size(C,1),size(steering_vector,2),size(rsp,1));
rsp = permute(rsp,[1 3 2]);
rsp = rsp(:,:,1:snapshots);
guard_bandwidth = floor(guard_bandwidth/2);
for ii = 1:size(rsp,1)

    % estimate spatial-doppler correlation matrix for each range using N snapshots
    % around given range value with a given guard band
    start_left = max([ii-(L+guard_bandwidth) 1 ]);
    end_left = max([ii - guard_bandwidth-1 1]);
    start_right = min([ii + guard_bandwidth+1 size(rsp,1)-1]);
    end_right = min([ii+(L+guard_bandwidth) size(rsp,1)-1]);
    
    % create interference plus noise covariance matrix using size a 3
    % weight doppler filter
    if(~max_power)
        s = (rsp([start_left:end_left start_right:end_right],:,:));
        s = permute(reshape(s(:),size(s,1),size(steering_vector,1),size(steering_vector_doppler,1),[]),[1 4 2 3]);
        s = reshape(s(:),[],size(steering_vector,1) * size(steering_vector_doppler,1));
        R = (s' * s)/size(s,1);
    else

         s = (rsp([start_left:end_left],:,:));
        s2 = abs(s(:,:,floor(size(s,3)/2)));
        [mv iii]=  sort(max(s2,[],2));
        s = s(iii(end:-1:max([end-points_kept 1])),:,:);
       % s = s(iii(end),:,:);
        s = permute(reshape(s(:),size(s,1),size(steering_vector,1),size(steering_vector_doppler,1),[]),[1 4 2 3]);
        s = reshape(s(:),[],size(steering_vector,1) * size(steering_vector_doppler,1));
        sl = s;
        s = (rsp([start_right:end_right],:,:));
        s2 = abs(s(:,:,floor(size(s,3)/2)));
        [mv iii]=  sort(max(s2,[],2));
        s = s(iii(end:-1:max([end-points_kept 1])),:,:);
        %s = s(iii(end),:,:);
        s = permute(reshape(s(:),size(s,1),size(steering_vector,1),size(steering_vector_doppler,1),[]),[1 4 2 3]);
        s = reshape(s(:),[],size(steering_vector,1) * size(steering_vector_doppler,1));
        sr = s;
        s = [sl;sr];
        R = (s' * s)/size(s,1);
    end
%         R = zeros(size(C,1),size(C,1));
%     %for jj = [1:size(steering_vector_doppler,1):15 24:size(steering_vector_doppler,1):snapshots-size(steering_vector_doppler,1)]
%    for jj = 1:size(steering_vector_dopplr,1):(snapshots)
% 
%         s = rsp([start_left:end_left start_right:end_right],:,jj:(jj-1)+size(steering_vector_doppler,1));
% 
%         s = reshape(s(:),size(s,1),[]);
%         R = R + (s' * s);
%    end

    % solve for weights
    x = linsolve(R+ sigma * mean(diag(R))  * eye(size(R)),C,opts);
    n =real(dot(C,x));
    w(:,:,ii) =   x./ sqrt(repmat(n,size(x,1),1));


end
%     R = zeros(size(C,1),size(C,1));
%     %for jj = [1:size(steering_vector_doppler,1):15 24:size(steering_vector_doppler,1):snapshots-size(steering_vector_doppler,1)]
%    for jj = 1:size(steering_vector_doppler,1):(snapshots)
% 
%         s = rsp([start_left:end_left start_right:end_right],:,jj:(jj-1)+size(steering_vector_doppler,1));
% 
%         s = reshape(s(:),size(s,1),[]);
%         R = R + (s' * s);
%    end
