for ii = size(track,3)+1
    %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
    track(:,:,ii) = track(:,:,2);
    track(2,:,ii) =[1:size(track,2)]*100;
    
    %track(1:3,:,ii+2) = 0;
    track(1,:,ii) =track(1,:,ii)+3000+12491.35;
    track(2,:,ii) =track(2,:,ii);
    track(3,:,ii) =track(3,:,ii);
    
end
reflextivity(end+1) = reflextivity(end)*5;