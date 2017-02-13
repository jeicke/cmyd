function [maxDelay, minDelay ] = finddelays( track,geometry ,C, pulselength)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
maxDelay = 0;
minDelay = 1E6;
count = 1;
if(~isempty(geometry))
    %for jj = 1:size(geometry,1)
        for ii = 1:size(track,2)
            
            
            xp = squeeze(track(1,ii,3:end)) - track(1,ii,2);
            yp = squeeze(track(2,ii,3:end)) - track(2,ii,2);
            zp = squeeze(track(3,ii,3:end)) - track(3,ii,2);
            range_multi_ant = sqrt(xp(:).^2 + yp(:).^2 + zp(:).^2);
            
            xp = squeeze(track(1,ii,3:end)) - track(1,ii,1);
            yp = squeeze(track(2,ii,3:end)) - track(2,ii,1);
            zp = squeeze(track(3,ii,3:end)) - track(3,ii,1);
            range_track_multi = (sqrt(xp(:).^2 + yp(:).^2 + zp(:).^2));
            
                maxDelay = max([maxDelay max(range_multi_ant(:) + range_track_multi(:))/C]);
                minDelay = min([minDelay min(range_multi_ant(:) + range_track_multi(:))/C]);
            
            delay(count,:) = (range_multi_ant(:) + range_track_multi(:))/C;
            count = count + 1;
        end
   % end
else
    
    for ii = 1:size(track,2)
        
        
        xp = squeeze(track(1,ii,3:end)- track(1,ii,2)) ;
        yp = squeeze(track(2,ii,3:end)- track(2,ii,2)) ;
        zp = squeeze(track(3,ii,3:end)- track(3,ii,2)) ;
        range_multi_ant = sqrt(xp(:).^2 + yp(:).^2 + zp(:).^2);
        
        xp = squeeze(track(1,ii,3:end)) - track(1,ii,1);
        yp = squeeze(track(2,ii,3:end)) - track(2,ii,1);
        zp = squeeze(track(3,ii,3:end)) - track(3,ii,1);
        range_track_multi = (sqrt(xp(:).^2 + yp(:).^2 + zp(:).^2));
        
            maxDelay = max([maxDelay max(range_multi_ant(:) + range_track_multi(:))/C]);
            minDelay = min([minDelay min(range_multi_ant(:) + range_track_multi(:))/C]);
        
        delay(count,:) = (range_multi_ant(:) + range_track_multi(:))/C;
        count = count + 1;
        
    end
end
x = delay;
minDelay = min(x(:));
maxDelay = max(x(:));
if(exist('pulselength','var')&&~isempty(pulselength))
    minDelay = minDelay;% - 1*(pulselength);
    maxDelay = maxDelay;% + 1*(pulselength);
end
end


