function SAresp = subarray(numel,spacing,steer)
%    numel
%    for spacing = .16:.02:.56
%        spacing
        posit = -((numel-1)/2)*spacing : spacing: (numel-1)/2*spacing;
        posit = posit';
        wts = chebwin(numel,40);
        wts = wts/sum(wts);
 %       for steer = 0:.04:1.04
 %           steer
             SAresp = calcBP(posit,wts,steer);
 %       end
 %       keyboard
 %   end
