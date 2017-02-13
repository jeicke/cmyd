function  [cpatime, cparange, velocity, angle, range, blockedTime, specularRange] = plottrack( track,step,movie,plotit)
%plottrack small test program to plot a track
if(isempty(step)||~exist('step','var'))
    step = 1:size(track,2);
end
if(isempty(plotit)||~exist('plotit','var'))
    plotit = true;
end
if(plotit)
    if(~movie)
        hold on
        plot(track(1,step ,1),track(2,step ,1),'rs',track(1,step ,2),track(2,step ,2),'go','markersize',16)
        for ii = 3:size(track,3)
            plot(track(1,step ,ii),track(2,step ,ii),'bd','markersize',3)
        end
       % legend('CELL TOWER','HANDSET','SCATTER POINT','location','west')
        grid on;
            xlabel('M')
            
            ylabel('M')
    else
        maxx = max(track(1,2:end));
        minx = min(track(1,2:end));
        
        maxy = max(track(2,2:end));
        miny = min(track(2,2:end));
        for jj = 1:size(track,2)
            for ii = 3:size(track,3)
                plot(track(1,step ,1),track(2,step ,1),'rs',track(1,step ,2),track(2,step ,2),'go','markersize',16)
                hold on;
                plot(track(1,jj ,ii),track(2,jj ,ii),'bd','markersize',8)
                
                
                
                
            end
            
            axis([minx maxx  miny maxy])
            hold off;
            grid on;
            xlabel('M')
            
            ylabel('M')
            drawnow;
        end
        
        
    end
end
r =  sqrt( (track(1,:,3)-track(1,:,2)).^2 + (track(2,:,3)-track(2,:,2)).^2 +(track(3,:,3)-track(3,:,2)).^2);
t = track(4,:,3);
tx = 0:(t(2)-t(1))/1000:t(end);
rx = interp1(t,r,tx);
[h i] = min(rx);
cpatime = tx(i);
cparange = rx(i);
r =  sqrt( (track(1,:,3)-track(1,:,1)).^2 + (track(2,:,3)-track(2,:,1)).^2 + (track(3,:,3)-track(3,:,1)).^2);
t = track(4,:,3);
tx = 0:(t(2)-t(1))/1000:t(end);
rx1 = interp1(t,r,tx);
[h i] = min(rx+rx1);
cpatime = [cpatime tx(i)];
cparange = [cparange h];
specularRange = rx(i);
blockedTime(1) = 1E9;
blockedTime(2) = -1;
for ii = 3:size(track,3)
    r =  sqrt( (track(1,:,ii)-track(1,:,2)).^2 + (track(2,:,ii)-track(2,:,2)).^2 +(track(3,:,ii)-track(3,:,2)).^2);
    t = track(4,:,3);
    tx = 0:(t(2)-t(1))/100:t(end);
    rx = interp1(t,r,tx);

    r =  sqrt( (track(1,:,ii)-track(1,:,1)).^2 + (track(2,:,ii)-track(2,:,1)).^2 + (track(3,:,ii)-track(3,:,1)).^2);
    t = track(4,:,3);
    tx = 0:(t(2)-t(1))/100:t(end);
    rx1 = interp1(t,r,tx);
    [h i] = min(rx+rx1);
    blockedTime(1) = min([blockedTime(1)  tx(i)]);
    blockedTime(2) = max([blockedTime(2)  tx(i)]);
  
end
[h i] = min(abs(t-cpatime(1)));
i1 = max([1 i-10]);
i2 = min([i+10 size(track,2)]);
dx = median(gradient((track(1,i1:i2,3))));
dy = median(gradient((track(2,i1:i2,3))));
dz = median(gradient((track(3,i1:i2,3))));
dt = (track(4,2,3)-track(4,1,3));
angle = atan2(dy,dx) * 180/pi-90;
if(angle<0)
    angle = angle + 360;
end
velocity = norm([dx dy dz]/dt);
range = r;
end
