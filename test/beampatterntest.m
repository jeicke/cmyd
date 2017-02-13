Fc = 9410E6;
C= 299792458;
arrayLength = 1.3;
lamda = C/Fc;
bw = 1.27 * lamda/arrayLength
width =  1/(25 * pi/180/(1.27 * lamda));

[geometry w h]= loadgeometry('rectangle',[],arrayLength,width,C/Fc);


[power az de parameters] = beampattern(Fc,[0],[-90:90],ones(length(geometry),1),[], [],geometry,[],C);