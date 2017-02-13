function [ c] = scatterers( m,k )
%scatterers calculates set of multipath scatterers
%   Detailed explanation goes here
for ii = 1:length(k)
    
    [c(ii,:) ] = findpoint(m,rand(1) * 360,norm(m) + k(ii) );
end
end

