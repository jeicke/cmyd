for ii = 1:size(p,2)
    ii
    %hold off;
   % plot((p(100:200,ii)))
   % hold on;
    plot((p2(940:1000,ii)),'r')
    axis([0 70 0 16])
    drawnow
    pause(.01)
end