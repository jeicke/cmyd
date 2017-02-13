%%
angle = [0:20:360]
count = 1;
clear psf B;
Ppsft = Ppsf(300:700,400:700,1);
for ii = 1:length(Pxx)-1
    psf(:,:,count) = gpuArray(double(abs(imrotate(Ppsft,angle(ii),'crop'))));
    B(:,:,count) = gpuArray(double(abs(Pxx{ii})));
    imagesc(B(:,:,count))
    drawnow

    count = count + 1;
end
g = B;
gn = B;
f = mean(B,3);
fn = mean(B,3);
figure
imagesc(f)
figure
for ii = 1:1000
    z = zeros(size(f));
    for jj = 1:5%size(psf,3)
        jj
        z = z + imfilter(g(:,:,jj)-gn(:,:,jj),psf(:,:,jj).^2/1000);
    end
    fn = fn+  z;
    imagesc(fn)

    drawnow
    for jj = 1:5%size(psf,3)
        gn(:,:,jj) = imfilter(fn,psf(:,:,jj)) ;
    end

end