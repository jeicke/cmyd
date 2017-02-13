rngspa= 299792458/(sampling_rate)/2;
ranges = [beam_transmit(5,1)*1000:rngspa:beam_transmit(4,1)*1000];
v = prf * c/operating_frequency/4 * ([-50:50]/50);
p = 20*log10(abs(squeeze(D(:,:,1))));
imagesc(v,ranges/1000,p);
colorbar
xlabel('RADIAL VELOCITY (m/s)');
ylabel('RANGE (KM)')


p = 20 * log10(abs(squeeze(max(Dbf2,[],2))));
%caxis([max(max(p))-100 max(max(p))])
figure;
count2 = 1;
count1 = 1;

for range = ranges
    for az = bearings
        
        X(count2,count1) = range * sind(az);
        Y(count2,count1) = range * cosd(az);
        count1 = count1 + 1;
    end
    count1 = 1;
    count2 = count2 + 1;
end
pcolor(X/1000,Y/1000,double(p));
shading flat;
axis equal
axis tight
xlabel('Y KM')
ylabel('X KM')
colorbar
caxis([max(max(p))-30 max(max(p))])

figure
p = 20 * log10(abs(squeeze(max(Dbf,[],2))));
pcolor(X/1000,Y/1000,double(p));
shading flat;
axis equal
axis tight
xlabel('Y KM')
ylabel('X KM')
colorbar
caxis([max(max(p))-23 max(max(p))])