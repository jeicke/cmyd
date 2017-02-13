%BEAMPATTERN computes a beampattern
% frequency           -- frequency at which to compute beampattern in Hz
%                        (scalar)
% az                  -- 1xN set of angles to calculate beampattern
% de                  -- 1xM set of angles to calculate beampattern
% w                   -- optional 1XSENSOR set of weights to use in computing
%                        beampattern
% parameters          -- optional structure containing stored parameters
%                        and computations to compute beampattern can
%                        contain...
%                        *parameters.Vs - speed of sound in units desired
%                         default 1524 m/s
%                        *parameters.SLx - SENSORX3 geometery of array in
%                         units compatable with parameters.Vs
%                        *parameters.theta_source if w is not defined will
%                         make a test source at the azimuthal angle specified in
%                         degrees
%                        *parameters.phi_source if w is not defined will
%                         make a test source at the D/E angle specified in
%                         degrees
% window               -- optional SENSORX1 window to apply to weights
% geometry             -- SENSORX3 matrix of sensor geometry in units
%                         compatable with Vs
% source               -- if w is not defined a 1x2 vector [AZ D/E] of
%                         source location in degrees
% Vs                   -- Velocity of sound in units compatable with
%                         geometery
% outputs
% power                -- MXN matrix of output power levels (linear scale)
% az
% de 
% parameters 
% ***below are optional commands
function [power, az, de, parameters] = beampattern(frequency,az,de,w,parameters,window,geometry,source,Vs,cr)

if(1)
    cast = @single;
else
    cast = @GPUsingle;
end
inverse_cast = @single;


usecr = true;
if ~exist('cr', 'var') || isempty(cr)
    usecr = false;
end

if(~exist('parameters')||isempty(parameters)||~isfield(parameters,'Vs'))
    if(~exist('Vs','var')||isempty(Vs))
       Vs = 1524; % 1524 is STAFAC speed of sound in m/s 
    end
else
    Vs = parameters.Vs;
end

if(~exist('parameters')||isempty(parameters)||~isfield(parameters,'SLx'))
    if(isempty(geometry)||~exist('geometry','var'))
        lengtha = 100;
        step = .02;
        geometry(:,1) = [-lengtha :step:lengtha ]';
        geometry(:,2) = zeros(size([-lengtha :step:lengtha ]'));
        geometry(:,3)= zeros(size([-lengtha :step:lengtha]'));
    end
    SLx = cast(geometry);
    parameters.SLx = SLx;
end
if(~exist('window','var')||isempty(window));
    window = cast(ones( size(parameters.SLx,1),1));
end
if(isfield(parameters,'subarray_geometry'))
    parameters_sub.Vs =Vs;
    parameters_sub.SLx = parameters.subarray_geometry;
    if(isfield(parameters,'subarray_theta_source'))
        parameters_sub.theta_source =parameters.subarray_theta_source;
        parameters_sub.phi_source=parameters.subarray_phi_source;
    else
        parameters_sub.theta_source = parameters.theta_source;
        parameters_sub.phi_source = parameters.phi_source;
    end
    if(isfield(parameters,'subarray_window'))
        subarray_window = parameters.subarray_window ;
    else
        subarray_window = [];
    end
    [theory_matrix] = beampattern(frequency,az,de,[],parameters_sub,subarray_window * length(subarray_window));
    
else
    theory_matrix = 1;
end
SLx = cast(parameters.SLx);
wavenumber = frequency/Vs;
sp_2pi = 2* pi* wavenumber ;

%[theta,phi] = meshgrid(az,de);
theta = (az(:)') * pi/180;
phi = (de(:)') * pi/180;
step_size = 100000;

%get index
start_index = 1;
end_index = step_size;
if(end_index>length(theta))
    end_index = length(theta);
end
index = start_index:end_index;
%kernal = zeros(length(geometry),length(index));
%step through look direction vector
power = cast(zeros(1,length(theta)));
sp_2pi  = cast(sp_2pi);
if(isfield(parameters,'usestored')&&parameters.usestored==true)
    if(~isfield(parameters,'savedE'))
       
        r = cast([sin(theta) .* cos(phi); cos(theta) .* cos(phi); sin(phi)]);
        parameters.savedE = exp( ( sp_2pi * SLx*r ) * 1i);
    end
    usestored = true;
else
    usestored = false;
end
if(~exist('w','var')||isempty(w))
    if(isfield(parameters,'theta_source')&&isfield(parameters,'phi_source'))
        
        w =  exp( ( sp_2pi * SLx* [sin(parameters.theta_source * pi/180) .* cos(parameters.phi_source* pi/180); cos(parameters.theta_source* pi/180) .* cos(parameters.phi_source* pi/180); sin(parameters.phi_source* pi/180)]) * 1i);
    else
        if(~exist('source','var')||isempty(source))
            
            w =  exp(  sp_2pi * SLx* cast([0; 1; 0])* 1i);
        else
            w =  exp( ( sp_2pi * SLx* [sin(source(1) * pi/180) .* cos(source(2)* pi/180); cos(source(1)* pi/180) .* cos(source(2)* pi/180); sin(source(2)* pi/180)]) * 1i);
        end
    end
    
    w = w.*window/sum(window.^2);
else
    w = single(w(:));
    w = cast(w) .* cast(window);
end


if(~(usestored))
    e = exp( ( sp_2pi * SLx*[sin(theta) .* cos(phi); cos(theta) .* cos(phi); sin(phi)] ) * 1i);
else
    e =  parameters.savedE(:,index);
end

% if usecr
%     Cr = crystal_response(theta, phi, cr)';
%     e = e.*sqrt(Cr);
% end

power =  abs(w'*e).^2;

power = inverse_cast(power);
power= (reshape(power,length(de),[]));
power = power.*theory_matrix;



