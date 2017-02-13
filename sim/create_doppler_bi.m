function [clutter_doppler, EL] = create_doppler_bi(AZ,velocity,range,frequency,altitude)




[EL graz]= sargmtiangles(altitude,sqrt(range.^2 + (altitude * 1000)^2)/1000);
EL= EL * pi/180;

A = sin(AZ);
B = cos(AZ);
C = cos(EL); 
D = sin(EL);

clutter_doppler = zeros(size(AZ,2),size(EL,2));

%2*target_velocity*frequency * 1E6/299792458
%clutterarea = sqrt(2 * range.^2 * (1 - cosd(resolution))) * rngspa ;
for ii = 1:length(D)
    clutter_doppler(:,ii) =   velocity' * [A.*C(ii);B.*C(ii);ones(size(A)) * D(ii)] * frequency/299792458;
end


