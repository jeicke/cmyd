clear track;
% receiver and transmitter locations
angle = 45;
receiverVelocity = 30* [sind(angle) cosd(angle) 0];
receiverPosition  =[cosd(angle) sind(angle) 0];
receiverPosition = receiverPosition-receiverVelocity * time/2;
receiverPosition(3) = 10E3;
transmitterPosition =[-20 0 20]*1000;

 transmitterVelocity =[0 0 0];
%  transmitterPosition = receiverPosition;
%  transmitterVelocity =receiverVelocity;
% make transmitter and receiver tracks
[ track] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );
targetVelocity = [10 0 0];

[tgttrack] =maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition );