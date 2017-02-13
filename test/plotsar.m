function [ output_args ] = plotsar( sar,x,y,dbdown )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
MMM = reshape(single(sar),length(y),[]);
p = 20*log10(abs(single((MMM))));
imagesc(x,[y(1) y(end)],(p))
set(gca,'YDir','normal')
caxis([max(p(:))-dbdown max(p(:))])
xlabel('X (m)')
ylabel('Y (m)')
colorbar;
 axis equal
 axis image
end

