function [elevation,grazing]=SARGMTIangles(height,slant_range)

rearth=6446*4/3; %4/3 earth approximation

raircraft=rearth+height;
cang=acos(-(slant_range.^2-raircraft^2-rearth^2)/(2*raircraft*rearth));
elevation=pi/2-asin(sin(cang)*rearth./slant_range);
%grazing=asin(sin(cang)*raircraft/slant_range)
grazing=pi-cang-(pi/2-elevation)-pi/2;
grazing=grazing*180/pi;
elevation=-elevation*180/pi;