% steering_vector Computes steering vector for a given geometry, bearing
% and wave number

% steering_vector(geometry, AZ, EL, K)
% is the steering vectors for the Nx3 geometry matrix at the given Azimuth
% and Elevation angles (in radians) for the K equal to 1/wave length
% steering_vector(geometry, AZ, EL, K,cal)
% cal is a vector Nx1 or NxN matrix of complex calibration values
function [V] = steering_vector(geometry, AZ, EL, K,cal)
if(false)
    [Beam_phase]=Beamsteer_Phase(AZ * 180/pi,true).';
    V = exp(i*Beam_phase);
else


    %V =exp((2 * pi * geometry * [sin(AZ+pi/2) .* cos(EL); cos(AZ+pi/2).* cos(EL); sin(EL)].*K)*i);
    V =exp((2 * pi * geometry * [sin(AZ) .* cos(EL); cos(AZ).* cos(EL); sin(EL)].*K)*-1i);
end
if(exist('cal','var')&&~isempty(cal)&&size(cal,2)==1)

    V=(V.*kron((cal),ones(1,size(V,2))));
elseif(exist('cal','var')&&~isempty(cal))
    V=cal*V;
end