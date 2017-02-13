clear all;
width = .2; %in m
height = .2; % in m
res = 1;

%% operating parameters
operating_frequency = 5E9;  % in Hz
c = 299792458; % speed of light in m/s

%% calculate beamsize in degrees
wavelength = c/operating_frequency;
beamfactor=17/20; %factor to convert to 3dB beamwidth

azbw = asind(2 * wavelength/width) * beamfactor;
elbw = asind(2 * wavelength/height) * beamfactor;
%% set if looking up or down
% if true we need to use grazing angles to steer beams, otherwise 
looking_down = true;


%% create example array 1/2 wavelength spaced
count = 1;
for w = [-width/2:wavelength/2:width/2]
    for h = [-height/2:wavelength/2:height/2]
        geometry(count,1) =w;
        geometry(count,2) = 0;
        geometry(count,3)= h;
        count = count +1;
    end
end
%%
theta = 1;
phi = 0;
kernal =  geometry * [sin(theta) * cos(phi); cos(theta) * cos(phi); sin(phi)] * -i * 2 * pi;
e = exp(kernal *  operating_frequency / c);

%%
nSens = size(geometry,1);
z = repmat(geometry,nSens,[]);

index = repmat([1:nSens],nSens,1);
index = index(:);
d = sqrt(sum( (z-geometry(index,:)).^2,2));
A = 4 * pi * sinc(2 * pi * operating_frequency/ c*d/pi);
A = reshape(A,nSens,nSens).';
S2 = .05;
%Ri = 

%%
B = inv(A +S2 * mean(diag(A)) * eye(size(A))+50000);
weights = B*e/(e'*B*e);
HGA.theta_source =0;
HGA.phi_source = 0;
HGA.SLx = geometry;
resolution = 1;
[theory_matrix az de HGA] = calcBeamPattern(operating_frequency,[0:resolution :360],[-90:resolution:89],weights,HGA);