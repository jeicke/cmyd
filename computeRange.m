function [r az de] = computeRange(track, times,start)
%compute range for all phones for times in vector time

az = [];
de = [];
%get range at each sample
%r = zeros(size(sourcePos,2),size(times,2));

x = interp1q(track(4,:,1)',-bsxfun(@plus,track(1,:,start),-squeeze(track(1,:,:))')',times');
y = interp1q(track(4,:,1)',-bsxfun(@plus,track(2,:,start),-squeeze(track(2,:,:))')',times');
z = interp1q(track(4,:,1)',-bsxfun(@plus,track(3,:,start),-squeeze(track(3,:,:))')',times');
r = sqrt(x.^2 + y.^2 + z.^2)';

az = atan2(x',y');
az(start,:) = [];

de = asin(z'./r);
de(start,:) = [];

r(start,:) = [];


% for ii = 1:size(sourcePos,2)
%
%
%     trackSource = bsxfun(@plus,track(1:3,:),sourcePos(:,ii));
%
%
%     x = interp1q(track(4,:)',trackSource(1,:)',times');
%     y = interp1q(track(4,:)',trackSource(2,:)',times');
%     z = interp1q(track(4,:)',trackSource(3,:)',times');
%
%     r(ii,:) = sqrt(x.^2 + y.^2 + z.^2)';
% end
% a = 1;
%

