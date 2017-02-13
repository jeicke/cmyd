function [DI PSG] = calculatePSG(resolution,BP,az,de,noise)
%[theory_matrix] = tbp(18,6,freq,true,[-180:1:180],[-90:90],[0
diacc = 0;
psgacc = 0;

for ii = 1:length(de)
   % for jj = 1:length(az)
        ringsum = sum(BP(ii,:),2);
        diacc  = diacc +  ringsum * cos(de(ii) * pi/180) * (resolution * pi/180) ^ 2;
        psgacc = psgacc + noise(ii) *  ringsum * cos(de(ii) * pi/180) * (resolution * pi/180) ^ 2;
   % end
end
%calculate DI and PSG
DI  = 10 * log10(4 * pi / diacc);
PSG = 10 * log10(4 * pi / psgacc); 

