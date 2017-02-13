%Platform geometry

% I'm assuming this velocity is platform velocity in m/s this equals 648
% kph
velocity=180;

% in km
height=5;

%elevation beamwidth
elbw=1/33; 

azbw=1/33;

minr=24;
maxr=100;
azcovbeams=70; 
beamfactor=17/20; %factor to convert to 3dB beamwidth
maxscanangle=sqrt(3)/2; %in sin angle space
scanlossexponent=2.4; %2 way scan loss exponent for cos(theta)^x model

%Radar Parameters
Gr=41.5;%+20*log10(3);
Gt=41.5;
NF=3;
losses=5;
Ptav=40;
lamd=0.03;
sigma=0;
SNR=25;
attn=2.8/300;
resolution=1;
resolution_area=10*log10(resolution^2);

%break up into elevation beams (Break points between beams)
[elev,graz,slant]=groundsurv(height,minr,maxr,elbw*beamfactor);
no_elev_beams=length(elev);
for k=(1:no_elev_beams-1)
   elbeampos(k)=asin(mean(sin(elev(k:k+1)*pi/180)))*180/pi; %elevation beam positions
end
slantbeampos=sargmtiangles2(height,elbeampos);

%Required dwell time for end of each beam position
Scalloploss=0; %Scalloping loss at far range
dwell=detectionrange2(Ptav,Gt,Gr,lamd,NF,sigma+resolution_area,losses+Scalloploss+attn.*slant,slant*1e3,SNR);

for k=1:no_elev_beams-1
   elevdwell(k)=sum(dwell(1:k));
   illuminated_range(k)=slant(k)-slant(k+1);
   sumillrange(k)=sum(illuminated_range(1:k));
end   
%scanlossfactor=1.25
%azdwell=elevdwell*azcovbeams*scanlossfactor;
azbeampos=-maxscanangle:(maxscanangle*2)/(azcovbeams-1):maxscanangle;
squintloss=cos(asin(azbeampos)).^scanlossexponent;
azimuthbeamdwell=kron(1./squintloss',dwell(1:no_elev_beams-1));
azdwell=sum(sum(azimuthbeamdwell));
azdwelllong=sum(azimuthbeamdwell(:,1)); %Uses tx sidelobes to fill in an near in beams
%azdwell=search time for energy limited configuration


%Required dwell time for SAR resolution
resdwell=sardwell(slant,velocity,lamd,resolution);
az=-0.866:1.732/(azcovbeams-1):0.866;
dialation=1./abs(cos(asin(az)))';
azresdwell=dialation*resdwell;
totresdwell=sum(sum(azresdwell));



   
