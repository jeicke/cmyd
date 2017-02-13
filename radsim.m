function Parameters = radsim(fileNameBase,time,track,snr,varargin)

%RADARSIM simulates a radar
% fileNameBase        -- base filename for simulation
% time                -- time in simulation in seconds
% track               -- tracks of transmitter, receiver and scatterers 
%                        as a 4xMxN matrix the first four  rows are x,y,z
%                        and time.  M is the number of entries in the track
%                        N is the tracks N=1 is the track of the
%                        transmitter.  N=2 is the track of the receiver
%                        N=3 and greater is the track of each scatterer.
%                        The simulation will trace the track from the
%                        transmitter to the scatterer and then to the
%                        receiver
% snr                 -- SNR of signals in dBw
%                        note the noise is not added you have to do it
%                        external to this program with the line
%                        1/sqrt(2) * (randn(size(data)) + 1i * randn(size(data)));
%                        (assuming the data is baseband)
% ***below are optional commands
% 'fc'                -- Carrier frequency 0 is non-baseband processing
%                        (default 10 GHz)
% 'fs'                -- Baseband sampling frequency (default 100 MHz)
% 'c'                 -- Speed of wave in m/s (default 299792458 m/s)
% 'timeseries'        -- Function handle to timeseries (default @tone)
% 'parameters'        -- Array of parameters to pass to timeseries function
%                        (default 5000 Hz)
% 'propogate'         -- function handle to propogate method
% 'clutter'           -- clutter mask matrix
% 'channels'          -- channels in simulation
% 'bidirectional'     -- true if transmission is radar like, (to target and
%                        back) otherwise, false (default false)
% 'rangecorrect'     -- if true, add in range propogation effects
% 'groundbounce'     -- if true, compute ground bounce
% 'supressdirect'    -- if true, supress direct path power
% 'sensorresponse'   -- structure that hold sensor response parameters
% 'geometry'          -- sensorx3 geometry of receive array
% 'transmitgeometry'  -- sensorx3 geometry of transmit array
% 'transmitorientation' -- orientation of transmitter  as (roll,pitch,yaw,t) matrix
% 'transmitshading'     -- sensorx3 matrix of shading for transmit array

% DEFAULT PARAMETERS
%carrier frequency in Hz
carrierFrequency = 0;
%sample rate at baseband
basebandSamplingRate = 2^18;
% speed of light in m/s
C = 299792458 ;

% power clipping in dB on receiver
powerClipping = 50;

nChannels = 1;
timeseries = @tone;
propogate = @propagatedatapulse2;
timeSeriesParameters = 9993.76;
rangeCorrect = false;
reflexi = [0 0 0];
sensorResponse = [];
groundBounce = false;
supressDirect = false;
geometry = [];
transmitGeometry = [];
transmitOrientation = [];
transmitShading = [];
% PARSE INPUT
fprintf('***Parsing Input...\n')
if(~ischar(fileNameBase))
    ME = MException('radarsim:fileNameBase', ...
        'File name base must be a character string');
    throw(ME);
end
%time of simulation in seconds
if(~exist('time','var')||isempty(time))
    time = 15;
else
    if(time <=0)
        ME = MException('radarsim:time', ...
            sprintf('Time must be greater than zero Time is %f',time));
        throw(ME);
    end
end

% transmitted snr in dB
if(~exist('snr','var')||isempty(snr))
    snr = 20;
end


% adjustable parameters
for ii = 1:length(varargin)
    if(ischar(varargin{ii}))
        switch lower(varargin{ii})
            case {'fc'}
                carrierFrequency = varargin{ii+1};
                if carrierFrequency < 0;
                    
                    ME = MException('radarsim:fc', ...
                        sprintf('fc must be >= 0'));
                    throw(ME);
                    
                end
            case {'propogate'}
                propogate = varargin{ii+1};
            case {'clutter'}
                clutter = [1 ; varargin{ii+1}];
           
            case 'fs'
                basebandSamplingRate = varargin{ii+1};
                if basebandSamplingRate < 0;
                    
                    ME = MException('radarsim:fs', ...
                        sprintf('fs must be > 0'));
                    throw(ME);
                    
                end
            case 'c'
                C = varargin{ii+1};
                if C < 0;
                    
                    ME = MException('radarsim:C', ...
                        sprintf('Speed of light must be  > 0'));
                    throw(ME);
                    
                end
            
            case 'rangecorrect'
                rangeCorrect = varargin{ii+1};
                if(rangeCorrect~=true&&rangeCorrect~=false)
                    
                    ME = MException('radarsim:range_correct', ...
                        sprintf('range_correct must be  true or false'));
                    throw(ME);
                    
                end
                
            case 'groundbounce'
                groundBounce = varargin{ii+1};
                if(groundbounce~=true&&groundbounce~=false)
                    
                    ME = MException('radarsim:ground_bounce', ...
                        sprintf('ground_bounce must be  true or false'));
                    throw(ME);
                    
                end
            case 'timeseries'
                
                timeseries = varargin{ii+1};
            case 'supressdirect'
                supressDirect = varargin{ii+1};
                if(supressDirect~=true&&supressDirect~=false)
                    
                    ME = MException('radarsim:supress_direct', ...
                        sprintf('supress_direct must be  true or false'));
                    throw(ME);
                    
                end
            case 'parameters'
                
                timeSeriesParameters = varargin{ii+1};
            case 'sensorresponse'
                sensorResponse = varargin{ii+1};
            case 'geometry'
                geometry = varargin{ii+1};
            case 'transmitgeometry'
                transmitGeometry = varargin{ii+1};
            case 'transmitorientation'
                transmitOrientation = varargin{ii+1};
            case 'transmitshading'
                transmitShading = varargin{ii+1};
            otherwise
                
        end
    end
end

if(isempty(geometry))
    geometry = [0 0 0];
end
if(isempty(transmitShading))
    transmitShading = ones(size(geometry));
end

% find system delays
[maxSystemDelay, minSystemDelay ] = finddelays( track,geometry ,C,timeSeriesParameters(2));
if(minSystemDelay<0)
    fprintf(1,'\nERROR minSystemDelay is less than zero increase range or decrease pulse length\n');
    Parameters  = [];
    return;
end
% remove zero power points
pos = reflexi';


if(rangeCorrect)

      theta = 2 * pi * rand(size(clutter));
  %  theta(1) = pi/2;
    clutterPower =clutter.*(sin(theta) + 1i * cos(theta));
else
   
    %clutterPower =   ones(1,length(clutter));
    %clutterPower = (10.^(snr(1:length(clutter))/20)) .*  ones(1,length(clutter));

    theta = 2 * pi * rand(size(clutter));
  %  theta(1) = pi/2;
    clutterPower =clutter.*(sin(theta) + 1i * cos(theta));
   % theta = 2 * pi * rand(size(clutterPower));
    %theta(1) = pi/2;
    %clutterPower =  clutterPower.*(sin(theta) + 1i * cos(theta));
end

if(supressDirect)
    clutterPower(1) = 0;
end
%clutterPower = clutterPower .* clutter';

% return parameters in structure
Parameters.timeseries = timeseries;
Parameters.track = track;
Parameters.clutterPower = clutterPower;
Parameters.basebandSamplingRate = basebandSamplingRate;
Parameters.carrierFrequency = carrierFrequency;
Parameters.C = C;
Parameters.time = time;
Parameters.fileNameBase = fileNameBase;
Parameters.maxSystemDelay = maxSystemDelay;
Parameters.minSystemDelay = minSystemDelay;
Parameters.timeSeriesParameters = timeSeriesParameters;
Parameters.receiveGeometry = geometry;
Parameters.rangeCorrect = rangeCorrect;
Parameters.sensorResponse = sensorResponse;
Parameters.transmitGeometry = transmitGeometry;
Parameters.transmitOrientation = transmitOrientation;

%delay data 
displayparameters(Parameters);

TSTART =tic;
fprintf('***Starting processing...\n')
blockSize = propogate(timeseries, track, clutterPower,basebandSamplingRate,...
    carrierFrequency,C ,time,fileNameBase,maxSystemDelay,...
    minSystemDelay,timeSeriesParameters,geometry,...
    rangeCorrect,sensorResponse,transmitGeometry,transmitOrientation,transmitShading);

Parameters.blockSize = blockSize;
save(sprintf('%s_%s',fileNameBase,'Parameters.mat'),'Parameters');
t = toc(TSTART);
displayparameters(Parameters);

if(t<1)
    fprintf(1,'Elapsed time is %9.3f miliseconds\n',t * 1000)
else
     fprintf(1,'Elapsed time is %9.3f seconds\n',t)
end

dt = time(end)-time(1);
fprintf(1,'***Megsamples Per Second %5.3f\n',((basebandSamplingRate * dt)/t)/1E6)

