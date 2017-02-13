function [slantrange]=SARGMTIangles2(height,elev)
%elev is negative if pointing at the earth
rearth=6446*4/3; %4/3 earth approximation

raircraft=rearth+height;
ang1=pi/2+elev*pi/180;
ang2=pi-asin(raircraft*sin(ang1)/rearth);
ang3=pi-(ang1+ang2);
slantrange=sin(ang3).*rearth./sin(ang1);