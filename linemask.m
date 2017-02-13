%% create ground return reflexivity mask
% load image to process
%length of pulse in seconds
layers = 1;
layerheight = [10 35]
clear sigma;
clear scatteredpoint;
count = 1;
for layer =1:layers
    
    reflexivity = ones(100,1);
    %reflexivity  = reflexivity  +1;
    pixelSpaceing = lamda/2; %each pixel size in m
    
    % calculate range spacing in meters
    
    
    x = [1:size(reflexivity ,1)];
    x = (x-mean(x)) * pixelSpaceing ;
    
    y = [1:size(reflexivity ,2)];
    y = (y-mean(y)) * pixelSpaceing ;
    
    z = ones(size(y)) *layerheight(layer);
    
    

    
    for ii = 1:size(reflexivity,1)
        for jj = 1:size(reflexivity,2)
            if(reflexivity(ii,jj)~=0)
                
                sigma(count) = (reflexivity(ii,jj));
                scatteredpoint(:,count) = [y(jj)+142 x(ii) z(jj)];
                %scatteredpoint(:,count) = [1000 70-count*30 z(jj)];
                count = count + 1;
            end
            
        end
    end
end
%scatteredpoint = [];
%scatteredpoint(:,1) = [1000 -80 0];


%sigma= [];
%sigma(1) = 1;

 
% f = 10GHz
% T = 15s
% prf  = 8500
% cbw = 180 
% fs =  400 mhz
% target 100 m wide
% 100 km range