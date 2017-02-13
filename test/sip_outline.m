%% create ground return reflexivity mask

%clear variables;
%close all;

layers = 1;
layerZs = [3 6 9 12 15 18];
invertImage = true;
dx = 2;

radarOffset = [0 6];

clear sigma;
clear spt;
count = 1;

for layer=1:layers
    
    image_filename = sprintf('ship%d.png',layer-1);
    image_filename = 'C:\Users\jeicke\Desktop\sprdn\ISAR code\bitmaps\test3.png';
    A = double(imread(image_filename)/255);
   % imagesc(A)
    %pause
    if(invertImage)
        A = 1-A;
    end
    
    A = mean(A,3);
    
    [m n] = size(A);
    
    pixelSpaceing = 2; %each pixel size in m
    
    % calculate range spacing in meters
    x = (1:m);
    x = (x-mean(x)) * dx ;
    
    y = (1:n);
    y = (y-mean(y)) * dx ;
    
    z = ones(size(y)) * layerZs(layer);
    
    % take into account radar offset %
    x = x - radarOffset(1);
    y = y - radarOffset(2);

    for ii=1:m
        for jj=1:n
            if(A(ii,jj) > 0)    
                sigma(count) = rand(size((A(ii,jj))));
                spt(:,count) = [y(jj) x(ii) z(jj)];
                count = count + 1;
            end
        end
    end
end
scatteredpoint = spt;
%% plot scatter points %%
figure;
plot3(spt(1,:), spt(2,:), spt(3,:), 'b.');
hold on;
plot3(0, 0, 16, 'ro');
hold off;
grid on;
axis equal;

