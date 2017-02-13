SAsteer = 37;
FAsteer = 57;
theta = -90:1:90;
spacing = 0.16: .01 : .56;
for i = 1:size(spacing,2)
    C(i,:) = fullarray(12,4,spacing(i),SAsteer,FAsteer);
end
figure
imagesc(theta,spacing,10*log10(C),[-80 0]);
%normax;
colorbar;
hold on;
%contour(theta,spacing,10*log10(C),[-0 -40],'k')
%contour(theta,spacing,10*log10(C),[-0 -40],'w:')
%ylabel('Spacing (wavelengths)');
%title('SA 12/4 SA steer 50 FA steer 57');