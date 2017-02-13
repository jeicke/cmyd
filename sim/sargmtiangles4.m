function [accessarea]=SARGMTIangles2(height,graz);
%elev is negative if pointing at the earth
rearth=6446*4/3; %4/3 earth approximation

raircraft=rearth+height;
ang1=pi/2+graz*pi/180;
ang2=asin(rearth*sin(ang1)./raircraft);
ang3=pi-(ang1+ang2);

h=rearth*(1-cos(ang3));
accessarea=2*pi*rearth*h;