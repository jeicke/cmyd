%% create ground return reflexivity mask
% load image to process
%length of pulse in seconds
invert_image = true;
reflexivity  = double(imread(image_filename ));
if(invert_image)
    reflexivity = ~reflexivity;
else
    reflexivity = (reflexivity(:,:,1)) *256;

end

%reflexivity  = reflexivity  +1;
pixelSpaceing = 10; %each pixel size in m

% calculate range spacing in meters
rngspa= 299792458/Fs;

x = [1:size(reflexivity ,1)];
x = (x-mean(x)) * pixelSpaceing ;

y = [1:size(reflexivity ,2)];
y = (y-mean(y)) * pixelSpaceing ;

z = ones(size(y)) *0;


count = 1;
clear sigma;
clear scatteredpoint;
for ii = 1:size(reflexivity,1)
    for jj = 1:size(reflexivity,2)
        if(reflexivity(ii,jj)~=0)
            
            sigma(count) = (reflexivity(ii,jj));
            scatteredpoint(:,count) = [x(ii)-50 y(jj) z(jj)];
            %scatteredpoint(:,count) = [1000 70-count*30 z(jj)];
            count = count + 1;
        end
        
    end
end
scatteredpoint = [];
scatteredpoint(:,1) = [1000 -80 0];


sigma= [];
sigma(1) = 1;

 
% f = 10GHz
% T = 15s
% prf  = 8500
% cbw = 180 
% fs =  400 mhz
% target 100 m wide
% 100 km range