function [slantrange,grange]=SARGMTIangles2(height,graz);
%elev is negative if pointing at the earth
rearth=6446*4/3; %4/3 earth approximation

raircraft=rearth+height;
ang1=pi/2+graz*pi/180;
ang2=asin(rearth*sin(ang1)./raircraft);
ang3=pi-(ang1+ang2);
slantrange=sin(ang3).*rearth./sin(ang2);
grange=ang3*rearth;