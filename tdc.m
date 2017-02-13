function [P subarray] =tdc( fileNameBase,Fs,Fc,C,x,y,z,coeff,method,time)
%Performs SAR time domain correlation
% fileNameBase        -- base filename for data to process
% 'fc'                -- Carrier frequency 0 is non-baseband processing
%                        (default 10 GHz)
% 'fs'                -- Baseband sampling frequency (default 100 MHz)
% 'C'                 -- Speed of wave in m/s (default 299792458 m/s)
% 'x'                 -- x component of processing grid (defaults to
%                        processed grid)
% 'y'                 -- y component of processing grid (defaults to
%                        processed grid)
% 'z'                 -- xz component of processing grid (defaults to
%                        processed grid)
% 'coeff'             -- coeff if true returns correlation coefficent
%                        rather than power
% 'method'            -- function pointer to interpolation method
% 'time'              -- [start end] time to process
cuda = false;
if(cuda)
    cast = @GPUsingle;
    inverse_cast = @double;
else
    cast = @double;
    inverse_cast = @double;
end
parallel = false;


P = [];
% get length of files
if(~exist('C','var')||isempty(C))
    C = 299792458;
end
if(~exist('coeff','var')||isempty(coeff))
    coeff = false;
end

[samples_received] = testread(sprintf('%s_%s',fileNameBase,'received.bin'),-1);
[samples_transmitted] = testread(sprintf('%s_%s',fileNameBase,'transmitted.bin'),-1);
if(samples_received~=samples_transmitted)
    fprintf(1,'TRANSMITTED AND RECEIVED SAMPLE NOT EQUAL\n');
end
track = load(sprintf('%s_%s',fileNameBase,'track.mat'));
if(~exist('x','var')||~exist('y','var'))
    x = track.xi;
    y = track.yi;
    z = track.zi;
end
if(exist('time','var')&&~isempty(time))
    start_sample = floor(time(1) * Fs);
    if(start_sample>samples_received)
        fprintf(2,'ERROR: TIME IS GREATER THAN TIME IN FILE');
        return;
        
    end
    end_sample = ceil(time(2) * Fs);
    if(end_sample>samples_received)
        end_sample = samples_received;
    end
else
    start_sample = 0;
    end_sample = samples_received;
end
x = cast(x);
y = cast(y);
z = cast(z);
blockLength =1000000;

totalBlocks = (end_sample-start_sample)/blockLength;

blocks = 0;
if(~exist('method','var')||isempty(method))
    interp_method =   @fastspline;%@fastslinear;%@fastspline;
else
    interp_method = method;
end

% calculate max system delay
max_system_delay = 0;
for ii = 1:size(track.radar_track,2)
    xp = x + track.radar_track(1,ii)/1000;
    yp = y + track.radar_track(2,ii)/1000;
    zp = z + track.radar_track(3,ii)/1000;
    rangesp = (sqrt(xp(:).^2 + yp(:).^2) + zp(:).^2) * 1000;
    xp = x + track.receiver_track(1,ii)/1000;
    yp = y + track.receiver_track(2,ii)/1000;
    zp = z + track.receiver_track(3,ii)/1000;
    range_rec = (sqrt(xp(:).^2 + yp(:).^2) + zp(:).^2) * 1000;
    max_system_delay = max([max_system_delay max(rangesp(:)+range_rec(:))/C]);
    
end


state = [];
state.maxSystemDelay = max_system_delay  *4;
P = cast(zeros(size(x,1),size(x,2)));

if(coeff)
    Co = cast(zeros(size(x,1),size(x,2)));
    s1 = cast(zeros(size(x,1),size(x,2)));
    s2 = cast(zeros(size(x,1),size(x,2)));
end
k = cast(-1i * 2 * pi  * Fc/C);
if( parallel)
    matlabpool open
end
data_received(1) = 1;
first = true;
max_system_delay  = floor(max_system_delay * Fs)/Fs;
while blocks<totalBlocks
    
    tStart = (blocks * blockLength + start_sample)/Fs;
    
    tEnd = (start_sample + blocks * blockLength + blockLength - 1)/Fs;
    
    if(tEnd * Fs>end_sample)
        tEnd =end_sample/Fs;
        blockLength = ceil((tEnd-tStart) * Fs);
    end
    
    if(first)
        if(tStart~=0)
            tOld =  tStart - max_system_delay-2/Fs;
            tOldEnd = tStart;
            offset = floor(tOld * Fs);

            
            [data_transmitted] = testread(sprintf('%s_%s',fileNameBase,'transmitted.bin'),...
                offset * 8, ceil(tOldEnd * Fs)-floor(tOld* Fs));
            times = (tOld:1/Fs:tOldEnd-1/Fs);
            state.create = true;
            l1 = min([length(times) length(data_transmitted)]);
            if(length(times)~= length(data_transmitted))
                fprintf(1,'WARNING VECTORS NOT EQUAL1\n');
            end
            times  = times(1:l1);
            [j, state] = signal( times,state,data_transmitted(1:l1).'  );
        end
    end
    first = false;
    offset = floor(tStart * Fs);

    tStart = offset/Fs;
    tEnd = (offset + blockLength-1)/Fs;
    [data_received] = testread(sprintf('%s_%s',fileNameBase,'received.bin'),...
        offset * 8,blockLength );
    [data_transmitted] = testread(sprintf('%s_%s',fileNameBase,'transmitted.bin'),...
        offset * 8,blockLength);
    
    times = (tStart:1/Fs:tEnd);
    state.create = true;
    l1 = min([length(times) length(data_transmitted)]);
    if(length(times)~= length(data_transmitted))
        fprintf(1,'WARNING VECTORS NOT EQUAL\n');
    end
    times  = times(1:l1);
    [j, state] = signal( times,state,data_transmitted(1:l1).'  );
    
    % interpolate track
    x_tr = cast(interp1(track.radar_track(4,:),track.radar_track(1,:),times,'linear','extrap'));
    
    y_tr = cast(interp1(track.radar_track(4,:),track.radar_track(2,:),times,'linear','extrap'));
    z_tr = cast(interp1(track.radar_track(4,:),track.radar_track(3,:),times,'linear','extrap'));
    
    x_r = cast(interp1(track.receiver_track(4,:),track.receiver_track(1,:),times,'linear','extrap'));
    
    y_r = cast(interp1(track.receiver_track(4,:),track.receiver_track(2,:),times,'linear','extrap'));
    z_r = cast(interp1(track.receiver_track(4,:),track.receiver_track(3,:),times,'linear','extrap'));
    interpolator.init = true;
    
    [v interpolator ] = interp_method(interpolator,Fs, state.times,state.savedTs,cast);
    
    times = cast(times);
    
    data_received = cast(data_received);
    subarray(:,blocks+1) = mean([x_r; y_r;z_r],2);
    %x_tr = x_tr-mean(x_tr);
   % x_r = x_r-mean(x_r);
    for ii =1:size(x,1)
        tic
        
        for jj =1:size(x,2)
            
            
            pos =[x(ii,jj); y(ii,jj); z(ii,jj)] * 1000;
            
            
           
            r_transmitter = sqrt( (x_tr+pos(1)).^2 + (y_tr+pos(2)).^2 + (z_tr+pos(3)).^2);
            r_receiver = sqrt( (x_r+pos(1)).^2 + (y_r+pos(2)).^2 + (z_r+pos(3)).^2);
            
            r = r_transmitter + r_receiver;
            
            time_interp = times- r/C;
            
            [timeSeries] = interp_method(interpolator,Fs, time_interp,[]);
            
            timeSeries = timeSeries .* exp(-k * r);
            % [timeSeries] = delayFilterNearField(state,track.radar_track,track.receiver_track,Fs,Fc,C , tStart,tEnd,pos/1000  ).';
            
            P(ii,jj) =  P(ii,jj) + timeSeries(1:length(data_received))*data_received;
            if(coeff)
                s1(ii,jj) = norm([s1(ii,jj) timeSeries(1:length(data_received))]);
                s2(ii,jj) = norm([s2(ii,jj) data_received.']);
                Co(ii,jj) =  abs(P(ii,jj))/sqrt(s1(ii,jj).^2*s2(ii,jj).^2);
            end
            
        end
        toc
        
    end
    blocks = blocks + 1;
    fprintf(1,'Percent done: %3.2f\n',100 * blocks/totalBlocks);
    p = 20*log10(abs(P));
    %toc
    %plot(p);
    if(coeff)
        Co(1,1)
    end
    if(size(x,1)~=1&&size(x,2)~=1)
        imagesc(inverse_cast(p))
        caxis([max(p(:))-100 max(p(:))]);
        drawnow;
    else
        plot(inverse_cast(p))
        grid;
        drawnow;
    end
    
end
if( parallel)
    matlabpool close
end
P = inverse_cast(P);
if(coeff)
    P = inverse_cast(Co);
end
save(sprintf('%s_%s',fileNameBase,'correlation.mat'),'P');
end

