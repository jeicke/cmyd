function [clutterarea ] = bistaticpatcharea(AZ,range,altitude,gamma)
clutter_ranges = range;
%find elevation and grazing angles using 4/3 earth approximation
dr = median(diff(clutter_ranges ));
slantRange = sqrt(clutter_ranges.^2 + (altitude * 1000)^2);
clutter_ranges1 = [clutter_ranges(1)-dr clutter_ranges(1:end-1)];
r1 = sqrt(clutter_ranges1.^2 + (altitude * 1000)^2);
[~, graz]= sargmtiangles(altitude,slantRange/1000);
%graz = 180/pi*atan2(altitude * 1000,clutter_ranges);
%graz = abs(graz);
%get angle difference
resolution = median(diff(AZ));
rngspa = median(diff(range));

clutterarea =  sqrt(2 * clutter_ranges.^2 * (1 - cosd(resolution))) * rngspa ./ cosd(graz) ;%range*rngspa./cos(graz*pi/180);


%clutterarea = resolution*(slantRange .*cosd(graz)-sqrt(slantRange .^2.*cosd(graz).^2-(slantRange .^2-r1.^2)));
clutterarea = clutterarea .* sind(graz).*(10.^(gamma/10));
                                      