% Computes the gradient backprojection for the super-resolution
% minimization function. This function implements the gradient of the
% level one norm between the projection of the estimated HR image and each
% of the LR images.
%
% Inputs:
% Xn - The current estimate of the HR image
% LR - A sequence of low resolution images
% Fmot - The tranlational motion for each LR frame
% Hpsf - The PSF function

%
% Outpus:
% The backprojection of the sign of the residual error
function G=gradientbackproject(Xn, LR, Hpsf)

% Note that shift and blur are comutative, so to improve runtime, we first
% filter the HR image


% Allocate shifted and decimated HR image
G = zeros(size(Xn));
if(0)
    [m,n] = size(Xn);
    
    % output size
    PQ(1)= 2*m- 1;
    PQ(2) = 2*n - 1;
    padC_m = ceil((m-1)./2);
    padC_n = ceil((n-1)./2);
    fftXn = fft2(Xn,PQ(1), PQ(2));
    for k=1:size(LR,3)
        
        % Shift and decimate HR image for each frame k
        HRsd = ifft2(fft2(Hpsf(:,:,k),PQ(1), PQ(2)).*fftXn);
        HRsd = HRsd(padC_m+1:m+padC_m, padC_n+1:n+padC_n);
        Gsign =  sign(HRsd-LR(:,:,k));
        pr = ifft2(fft2(flipud(fliplr(Hpsf(:,:,k))),PQ(1),PQ(2)).*fft2(Gsign,PQ(1),PQ(2)));
        
        G(:,:,k)=pr(padC_m+1:m+padC_m, padC_n+1:n+padC_n);
        
    end
else
    for k=1:size(LR,3)
        
        % Shift and decimate HR image for each frame k
        
        Gsign =  sign(imfilter(Xn,Hpsf(:,:,k),0)-LR(:,:,k));
        pr = imfilter(Gsign, flipud(fliplr(Hpsf(:,:,k))), 0);
        
        G(:,:,k)=pr;
        
    end
end
% Compute the sum over k of the backprojected gradient
G = sum(G, 3);