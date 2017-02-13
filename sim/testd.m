velocity =100;
height = 10;
clutter_range = [sqrt(2 * rearth * height + height^2)]
angle = 60;
theta = 0;
prf = 1500;
range = 70* 1000;
operating_frequency = 5E9

ur1 = c/(2 * (prf * 1))/1000
clear dr1
format bank
for ii = 1:ceil(clutter_range/ur1)
    frange = range * ii
    [phi ]= sargmtiangles(height,sqrt(frange .^2 + (height * 1000)^2)/1000);
   
    dr1(ii) = 2 * velocity * (sind(angle) * sind(theta) * cosd(phi) + cosd(angle) * cosd(theta) * cosd(phi))* operating_frequency/299792458;
end


ur2 = c/(2 * (prf * 2))/1000
clear dr2
if(range>(ur2 * 1000))
    range = range-(ur2 * 1000);
    
end
for ii = 1:ceil(clutter_range/ur2)
    frange = range * ii
    [phi ]= sargmtiangles(height,sqrt(frange .^2 + (height * 1000)^2)/1000)
    
    dr2(ii) = 2 * velocity * (sind(angle) * sind(theta) * cosd(phi) + cosd(angle) * cosd(theta) * cosd(phi))* operating_frequency/299792458;
end
