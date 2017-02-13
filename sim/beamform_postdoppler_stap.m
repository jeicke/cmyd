function [dbf w sinr] = beamform_postdoppler_stap(rsp,steering_vector,snapshots,guard_bandwidth,sigma,points_kept)
count = 1;
sinr = [];

% create kronicker product steering vector

L = ceil(snapshots/2);

%rsp = reshape(rsp(:),1,[]);
if(~exist('max_power','var'))
    max_power = false;
end
opts.POSDEF = true; opts.SYM = true;


rsp = permute(rsp,[1 3 2]);
guard_bandwidth = floor(guard_bandwidth/2);
for jj = 1:size(rsp,3)
    for ii = 1:size(rsp,1)
        
        % estimate spatial-doppler correlation matrix for each range using N snapshots
        % around given range value with a given guard band
        start_left = max([ii-(L+guard_bandwidth) 1 ]);
        end_left = max([ii - guard_bandwidth-1 1]);
        start_right = min([ii + guard_bandwidth+1 size(rsp,1)-1]);
        end_right = min([ii+(L+guard_bandwidth) size(rsp,1)-1]);
        
        % create interference plus noise covariance matrix using size a 3
        % weight doppler filter
       
            
            s = (rsp([start_left:end_left],:,jj));
            s2 = abs(s(:,floor(size(s,2)/2)));
            [mv iii]=  sort(max(s2,[],2));
            s = s(iii(end:-1:max([end-points_kept 1])),:);
            % s = s(iii(end),:,:);
            
            sl = s;
            s = (rsp([start_right:end_right],:,jj));
            s2 = abs(s(:,floor(size(s,2)/2)));
            [mv iii]=  sort(max(s2,[],2));
            s = s(iii(end:-1:max([end-points_kept 1])),:);
            %s = s(iii(end),:,:);
           
            sr = s;
            s = [sl;sr];
            R = (s' * s)/size(s,1);
      
        %         R = zeros(size(C,1),size(C,1));
        %     %for jj = [1:1:15 24:1:snapshots-1]
        %    for jj = 1:size(steering_vector_dopplr,1):(snapshots)
        %
        %         s = rsp([start_left:end_left start_right:end_right],:,jj:(jj-1)+1);
        %
        %         s = reshape(s(:),size(s,1),[]);
        %         R = R + (s' * s);
        %    end
        
        % solve for weights
        x = linsolve(R+ sigma * mean(diag(R))  * eye(size(R)),steering_vector,opts);
        n = real(dot(steering_vector,x));
        wt =   x./ sqrt(repmat(n,size(x,1),1));
        X = rsp(ii,:,jj).';
        dbf(ii,jj) = (wt' * X);
        sinr(ii,jj) = abs(wt'*steering_vector).^2/abs(wt' * R * wt);
        w(ii,jj,:)  = wt;
        
        
    end
end