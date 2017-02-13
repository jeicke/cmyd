function [P] =tdc_pulse( fileNameBase,prf,pulse_replica,track,geometry,x,y,z,method)
%Performs SAR time domain correlation
% fileNameBase        -- base filename for data to process=
% 'prf'               -- Pulse repetition frequency in Hz
% 'track'             -- Track of simulation objects
% 'x'                 -- x component of processing grid
% 'y'                 -- y component of processing grid
% 'z'                 -- z component of processing grid
% 'method'            -- function pointer to interpolation method
cuda = false;
if(cuda)
    cast = @GPUsingle;
    inverse_cast = @double;
else
    cast = @double;
    inverse_cast = @double;
end
parallel = false;

% get length of files

C = 299792458;

% read file header and get info on run
% Fs sample rate in Hz
% Fc carrier frequency in Hz
[samples_received, Fs, Fc, channels] = testread(sprintf('%s',fileNameBase),-1);
% find min and max delay of saved data
pulseLength = length(pulse_replica)/Fs;
[systemDelay(2), systemDelay(1)] = finddelays( track,geometry ,C,pulseLength);
%find samples in each pulse of processed data (based on how simulator
%works)
L = floor((systemDelay(2)-systemDelay(1))*Fs)+floor(length(pulse_replica)/2);
if(L<2*length(pulse_replica))
    L = 2*length(pulse_replica );
end
L = 2.^(nextpow2(L));
%block length is number of stored sample per pulse
blockLength =L-floor(length(pulse_replica)/2);
%calculate total number of pulses in file
totalBlocks = floor(samples_received/L)-1;
%
blocks = 0;
state = [];

%weighting  =  hannwindow(length(pulse_replica)) .';
% weight transmitted pulse with low sidelobe weighting and prepair it for
% interpolation

weighting  = chebwin(length(pulse_replica),70);
pulse_replica =  double(pulse_replica./norm(pulse_replica));
filterz  = (pulse_replica).* weighting;
if(~exist('method','var')||isempty(method))
    interp_method =   @fastspline;%@fastslinear;%@fastspline;
else
    interp_method = method;
end
interpolator.init = true;

s = zeros(L,1);
s(1:length(filterz))=filterz;
times = [0:(length(s)-1)]/Fs+systemDelay(1);
[v, interpolator ] = interp_method(interpolator,Fs,times,s.',cast);


% calculate carrier phase
k = cast(-1i * 2 * pi  * Fc/C);

first = true;
P = zeros(size(x));
blockTimeStep = round(Fs*1/prf)/Fs;

while blocks<totalBlocks
    
    offset = blockLength*blocks+1;
    
    tStart = offset/Fs;
    tEnd = (offset + blockLength-1)/Fs;
    
    [data_received] = testread(sprintf('%s',fileNameBase), offset * 8,blockLength );
    data_received = cast(data_received);
    %  data_received = data_received  + 1/sqrt(2) * (randn(size(data_received ))+ 1i * randn(size(data_received )));
    
    state.create = true;
    D = data_received;
    % interpolate track
    t = track(4,:,1);
    blockTime= blocks*blockTimeStep+(systemDelay(2)-systemDelay(1))/2 +systemDelay(1) ;
    %interpolate track to time for pulse
    x_tr = cast(interp1(t ,track(1,:,1),blockTime,'linear','extrap'));
    y_tr = cast(interp1(t ,track(2,:,1),blockTime,'linear','extrap'));
    z_tr = cast(interp1(t ,track(3,:,1),blockTime,'linear','extrap'));
    
    x_r = cast(interp1(t ,track(1,:,2),blockTime,'linear','extrap'));
    y_r = cast(interp1(t ,track(2,:,2),blockTime,'linear','extrap'));
    z_r = cast(interp1(t ,track(3,:,2),blockTime,'linear','extrap'));
    for ii =1:size(x,1)
        tic
        for jj =1:size(x,2)
            pos =[x(ii,jj); y(ii,jj); z(ii,jj)];
            r_transmitter = sqrt( (x_tr+pos(1)).^2 + (y_tr+pos(2)).^2 + (z_tr+pos(3)).^2);
            r_receiver = sqrt( (x_r+pos(1)).^2 + (y_r+pos(2)).^2 + (z_r+pos(3)).^2);
            r = r_transmitter + r_receiver;
            time_interp = times(1:blockLength)+r/C-systemDelay(1);
            [timeSeries] = interp_method(interpolator,Fs, time_interp,[]);
            
            timeSeries = timeSeries .* exp(-k * r);
            % [timeSeries] = delayFilterNearField(state,track.radar_track,track.receiver_track,Fs,Fc,C , tStart,tEnd,pos/1000  ).';
            
            P(ii,jj) =  P(ii,jj) + timeSeries*data_received.';
        end
    end
    
    blocks = blocks + 1;
    
    fprintf(1,'Percent done: %3.2f\n',min([100 100 * blocks/totalBlocks]));
    % p = 20*log10(abs(P));
    first = false;
end



