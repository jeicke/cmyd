for ii = 1:1000
    z = zeros(size(f));
    for jj = 1:size(psf,3)
        pf = fft2(g(:,:,jj)-gn(:,:,jj),PQ(1), PQ(2));
        pr = real(ifft2(pf.*psf(:,:,jj)));
        pr = pr(2:size(Ppsf{1},1)+1,2:size(Ppsf{1},2)+1);
        z = z +pr;
    end
    fn = fn+  z;
    imagesc(20*log10(abs(fn)))
    pause
    drawnow
    for jj = 1:size(psf,3)
        rx = real(ifft2(fft2(fn,PQ(1), PQ(2)).*psf(:,:,jj)));
        gn(:,:,jj) =  rx(2:size(Ppsf{1},1)+1,2:size(Ppsf{1},2)+1);
    end

end

