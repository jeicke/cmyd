function [ track] =maketrack(time,step,transmitterVelocity,receiverVelocity,transmitterStart,receiverStart)
%MAKETRACK makes a track
%   make track makes the tracks used in the simulation
% time                -- total time in simulation in units compatable with
%                        velocity
% step                -- time step between track points
% transmitterVelocity -- 1x3 column vector of transmitter platform velocity
% receiverVelocity    -- 1x3 column vector of receiver platform velocity
% transmitterStart    -- Start location of transmitter (default [0 0 0])
% receiverStart       -- Start location of receiver (default [0 0 0])
if(~exist('transmitterVelocity','var')||isempty(transmitterVelocity))
    transmitterVelocity = [0 0 0];
else
    
    if(size(transmitterVelocity,1)~=1&&size(transmitterVelocity,2)~=2)
        ME = MException('radarsim:radar_velocity', ...
            sprintf('velocity must be a 1x3 column vector'));
        throw(ME);
    end
end
%receiver velocity in m/s
if(~exist('receiverVelocity','var')||isempty(receiverVelocity))
    receiverVelocity = [0 0 0];
else
    
    if(size(receiverVelocity,1)~=1&&size(receiverVelocity,2)~=2)
        ME = MException('radarsim:receiver_velocity', ...
            sprintf('velocity must be a 1x3 column vector'));
        throw(ME);
    end
end

tpps = 1/step;


s = ([1: floor((time+2)*tpps)]-1)/tpps;
track(1,:,1) = s * transmitterVelocity(1) + transmitterStart(1);% - radar_velocity(1)/2  * time;%s([1:Npulses]-1)/prf] * (radar_velocity) + radar_position;
track(2,:,1) = s * transmitterVelocity(2) + transmitterStart(2);% - radar_velocity(2)/2  * time;
track(3,:,1) = s * transmitterVelocity(3) + transmitterStart(3);% - radar_velocity(3)/2  * time;
track(4,:,1) = s;

track(1,:,2) = s * receiverVelocity(1) + receiverStart(1);%s([1:Npulses]-1)/prf] * (radar_velocity) + radar_position;
track(2,:,2) = s * receiverVelocity(2) + receiverStart(2);
track(3,:,2) = s * receiverVelocity(3) + receiverStart(3);
track(4,:,2) = s;

end

