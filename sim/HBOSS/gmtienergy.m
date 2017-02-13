%Platform geometry
height=20;
elbw=1/10000; 
minr=25;
maxr=300;
snr=16.5;
ang_cov=120; %angular coverage in degrees
scanloss=20*log10(mean(cosd(0:(ang_cov/2)).^1.5)); %scanloss=0;

%break up into elevation beams (Break points between beams)
[elev,graz,slant]=groundsurv(height,minr,maxr,elbw);
no_elev_beams=length(elev);

%Compute the solid angle for each elevation beam position
h=cosd(elev)*elbw;

%Transmitter power for 1m^2 aperture for each elevantion beam position
Pt=radarpower(0,1,snr,3,0,5,10,slant*1e3);

%Scale transmitter power for solid angles and angular coverage
Pte=Pt+10*log10(h)+10*log10(ang_cov/360)-scanloss;

%undB and sum erergy
energy=sum(10.^(Pte/10))
