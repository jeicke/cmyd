function [theory_matrix az de HGA] = calcBeamPattern(frequency,ang_az,ang_de,w,HGA,window,mask)

theta = single(zeros(1,length(ang_az) * length(ang_de)));
phi = theta;
count = 1;
if(~exist('HGA')||isempty(HGA)||~isfield(HGA,'Vs'))
	Vs = 299792458; % in m/s 
else
    Vs = HGA.Vs;
end

if(~exist('HGA')||isempty(HGA)||~isfield(HGA,'SLx'))
    lengtha = 4.5/2;
    step = .01;
    geometry(:,1) = [-lengtha :step:lengtha ]';
     geometry(:,2) = zeros(size([-lengtha :step:lengtha ]'));
      geometry(:,3)= zeros(size([-lengtha :step:lengtha]'));
    SLx = single(geometry );size
    HGA.SLx = SLx;
end
if(~exist('mask','var')||isempty(mask));
    mask = ones( size(HGA.SLx,1),1);
end
if(isfield(HGA,'subarray_geometry'))
    HGA_sub.Vs =Vs;
    HGA_sub.SLx = HGA.subarray_geometry;
    if(isfield(HGA,'subarray_theta_source'))
        HGA_sub.theta_source =HGA.subarray_theta_source;
        HGA_sub.phi_source=HGA.subarray_phi_source;
    else
         HGA_sub.theta_source = HGA.theta_source;
         HGA_sub.phi_source = HGA.phi_source;
    end
    if(isfield(HGA,'subarray_window'))
        subarray_window = HGA.subarray_window ;
    else
        subarray_window = [];
    end
    [theory_matrix] = calcBeamPattern(frequency,ang_az,ang_de,[],HGA_sub,subarray_window * length(subarray_window));

else
    theory_matrix = 1;
end
SLx = HGA.SLx;
wavenumber = frequency/Vs; 
sp_2pi = 2* pi* wavenumber ;
for the  = 1:length(ang_az)
    for ph = 1:length(ang_de)
        theta(count) = single(ang_az(the));
        phi(count)  = single(ang_de(ph));
        count = count + 1;
    end
end
theta = theta * pi/180;
phi = phi * pi/180;
step_size = 1000;

%get index
start_index = 1;
end_index = step_size;
if(end_index>length(theta))
    end_index = length(theta);
end
index = start_index:end_index;
%kernal = zeros(length(geometry),length(index));
%step through look direction vector
power = single(length(theta));
sp_2pi  = single(sp_2pi);
if(isfield(HGA,'usestored')&&HGA.usestored==true)
    if(~isfield(HGA,'savedE'))
        HGA.savedE = exp( ( sp_2pi * SLx* [sin(theta) .* cos(phi); cos(theta) .* cos(phi); sin(phi)]) * 1i);
    end
    usestored = true;
else
    usestored = false;
end    
if(~exist('w','var')||isempty(w))
    if(isfield(HGA,'theta_source')&&isfield(HGA,'phi_source'))
        
        w =  exp( ( sp_2pi * SLx* [sin(HGA.theta_source * pi/180) .* cos(HGA.phi_source* pi/180); cos(HGA.theta_source* pi/180) .* cos(HGA.phi_source* pi/180); sin(HGA.phi_source* pi/180)]) * i);
    else

        w =  exp(  sp_2pi * SLx* [sin(0) .* cos(0); cos(0) .* cos(0); sin(0)] * i);
    end
    w = w.*mask/(size(SLx,1));
else
    w = single(w .* mask);
end
if(~exist('window','var')||isempty(window))
    while (start_index<=length(theta))
        if(~(usestored))
            e = exp( ( sp_2pi * SLx* [sin(theta(index)) .* cos(phi(index)); cos(theta(index)) .* cos(phi(index)); sin(phi(index))]) * 1i);
        else
            e =  HGA.savedE(:,index);
        end

        power(index) =  abs(w'*e).^2;
        % update index
        start_index = start_index + step_size;
        end_index = end_index + step_size;
        index = start_index:end_index;
        if(end_index>length(theta))
            index = start_index:length(theta);
        end


    end
else

    while (start_index<=length(theta))
        if(~(usestored))
            e = exp( ( sp_2pi * SLx* [sin(theta(index)) .* cos(phi(index)); cos(theta(index)) .* cos(phi(index)); sin(phi(index))]) * 1i);
        else
            e =  HGA.savedE(:,index);
        end

        power(index) =  abs((w.*window)'*e).^2;
        % update index
        start_index = start_index + step_size;
        end_index = end_index + step_size;
        index = start_index:end_index;
        if(end_index>length(theta))
            index = start_index:length(theta);
        end


    end
end
power= (reshape(power,length(ang_de),[]));
power = power.*theory_matrix;
az = ang_az;
de = ang_de;
theory_matrix  = power;

