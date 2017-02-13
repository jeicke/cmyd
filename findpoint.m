function [p ] = findpoint(m,theta,k )
%FINDPOINT Finds multipath delay scattering points that have a given
%arrival angle and time
% m             --1x3 location of source in m
% theta         --arrival angle at receiver
% k             --total path length in m

mn = norm(m);
phi = atan2(m(2),m(1))-pi/180 * theta;
r2 = (mn^2-k^2)/(2 * (mn * cos(phi)-k));
p(1) = cos(theta * pi/180) * r2;
p(2) = sin(theta * pi/180) * r2;
p(3) = 0;
r1 = norm(m-p);
if( abs(norm((r1 + r2)-k))>.1)
    p = [0 0 0];
    fprintf(1,'***warning no solution found\n');
end
end

