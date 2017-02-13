%%
angle = [0:20:360]
count = 1;
PQ = paddedsize(size(Ppsf));
clear psf psfg;
[m,n] = size(Ppsf);
[mb,nb] = size(Ppsf); 
% output size 
mm = m + mb - 1;
nn = n + nb - 1;
padC_m = ceil((mb-1)./2);
padC_n = ceil((nb-1)./2);
for ii = 1:length(Pxx)-1
    x1 = imrotate(Ppsf,360-angle(ii),'crop');

    psf(:,:,count) = fft2(abs(x1),mm,nn);
    psfg(:,:,count) = fft2(abs(x1).^2,mm,nn);
    
    psfg(:,:,count) =psf(:,:,count)/(max(max(abs(psf(:,:,count))))*5);
    psf(:,:,count) =psf(:,:,count)/(max(max(abs(psf(:,:,count)))));
    B(:,:,count) = abs(Pxx{ii}) ;
    figure(1)
    imagesc(x1)
    figure(2)
    imagesc( B(:,:,count))
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
%%
for ii = 1:1000
    z = zeros(size(f));
    for jj = 1:length(Pxx)-1
        pf = fft2(g(:,:,jj)-gn(:,:,jj),mm,nn);
        pr = ifft2(pf.*psfg(:,:,jj));
        pr = pr(padC_m+1:m+padC_m, padC_n+1:n+padC_n);
        z = z +pr;
        jj
    end
    fn = fn+  z;
    imagesc(fn);
    drawnow
    for jj = 1:length(Pxx)-1
        rx = ifft2(fft2(fn,mm,nn).*psf(:,:,jj));
        gn(:,:,jj) =  rx(padC_m+1:m+padC_m, padC_n+1:n+padC_n);
    end

end