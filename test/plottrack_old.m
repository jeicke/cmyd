function  plottrack( track,step)
%plottrack small test program to plot a track
if(isempty(step)||~exist('step','var'))
    step = 1:size(track,2);
end
plot(track(1,step ,1),track(2,step ,1),'rs',track(1,step ,2),track(2,step ,2),'go','markersize',16)
hold on
for ii = 3:size(track,3)
    plot(track(1,step ,ii),track(2,step ,ii),'bd','markersize',8)
end
legend('TRANSMITTER','RECIEVER','SCATTER POINT')
grid on;
xlabel('METERS')

ylabel('METERS')
end