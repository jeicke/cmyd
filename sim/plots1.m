%% plot_track
x = track.xyz(:,1);
y = track.xyz(:,2);
z = track.xyz(:,3);
figure
plot3(x,y,z);
hold on; grid on;
rotate3D

%% plot mic array
figure
hold on; box on;
for j=1:length(p)
    plot3(p(1,j), p(2,j), p(3,j), '.')
    text(p(1,j), p(2,j), p(3,j), num2str(j));
end

h = line([p(1,1) p(1,2)], [p(2,1) p(2,2)], [p(3,1) p(3,2)]);
h = line([p(1,2) p(1,3)], [p(2,2) p(2,3)], [p(3,2) p(3,3)]);
h = line([p(1,3) p(1,1)], [p(2,3) p(2,1)], [p(3,3) p(3,1)]);
h = line([p(1,4) p(1,1)], [p(2,4) p(2,1)], [p(3,4) p(3,1)]);
h = line([p(1,4) p(1,2)], [p(2,4) p(2,2)], [p(3,4) p(3,2)]);
h = line([p(1,4) p(1,3)], [p(2,4) p(2,3)], [p(3,4) p(3,3)]);
grid on
rotate3d

%% plot data
n = size(data_matrix,1);
figure
plot(data_matrix)
legend('Black','Blue','Red','Yellow')