pulse_replica=m_chirp(Parameters.timeSeriesParameters(1),Parameters.timeSeriesParameters(2),0,Parameters.basebandSamplingRate);
[h,i] = xcorr(pulse_replica);
h = abs(h);
[~,hc] = max(abs(h));

I = zeros(100,100);
I(50:60,50:60) = 1;
I(5,5) = 10;
I(60,60) = 1;
h = h/norm(h);
s = h(hc-length(I)/2:hc+length(I)/2-1);
psf = kron(s,s');

B = imfilter(I,psf) ;
g = B;
gn = B;
f = B;
fn = B;
%%
for ii = 1:1000
    fn = fn+ imfilter(g-gn,psf.^2/100) ;
    imagesc(fn)
    drawnow
    gn = imfilter(fn,psf) ;

end
%%
pulse_replica=m_chirp(200E6,Parameters.timeSeriesParameters(2),0,Parameters.basebandSamplingRate);
 image_filename = sprintf('wife.png',layer-1);
    
    A = double(imread(image_filename));
    A = A(1:48,1:48);
    I = A;
    [h,i] = xcorr(pulse_replica,pulse_replica.*chebwin(length(pulse_replica)));
h = abs(h);
[~,hc] = max(abs(h));
h = h/norm(h);
s = h(hc-length(I)/2:hc+length(I)/2-1);
s = zeros(size(s));
s(length(I)/2-1) = 1;
s(length(I)/2) = 1;
s(length(I)/2+1) = 1;
psf = kron(s,s');

B = imfilter(I,psf) ;
g = B;
gn = B;
f = B;
fn = B;
figure
imagesc(I)
figure
%%
for ii = 1:1000
     gn = imfilter(fn,psf) ;
    fn = fn+ imfilter(g-gn,psf.^2/100) ;
    imagesc(fn)
    drawnow
    error = norm(g-gn)

end
%%
addpath('C:\Users\jeicke\Desktop\sprdn\ISAR code\')
addpath('C:\Users\jeicke\Desktop\sprdn\ISAR code\test')
addpath('C:\Users\jeicke\Desktop\sprdn\ISAR code\bitmaps')
%%
pulse_replica=m_chirp(Parameters.timeSeriesParameters(1)*.5,Parameters.timeSeriesParameters(2),0,Parameters.basebandSamplingRate);
[h,i] = xcorr(pulse_replica);
h = abs(h);
[~,hc] = max(abs(h));
image_filename = sprintf('test4.png');
 A = double(imread(image_filename));
A =  xor(sum(abs(A),3),255)
A = A>0;
A = padarray(A,[100,100]);
I = double(A);
I = I(1:220,1:220);
figure;

imagesc(I)
 xlim([80 150])
    ylim([80 150])
h = h/norm(h);
s = h(hc-length(I)/2:hc+length(I)/2-1);
[H,i1] = max(s);
mask = zeros(size(s));
mask(i1-50:i1+50) = 1;

psfb = kron(s,s');
maskpsf= kron(mask,mask');
clear psf B;
count = 1;
figure
for angle = [0:5:90]
    p= imrotate(psfb,angle,'crop');
    psf(:,:,count) =   p/sum(abs(p(:)));
    B(:,:,count) = imfilter(I,psf(:,:,count),0) ;
    psf(:,:,count) = psf(:,:,count).*maskpsf;
    imagesc(B(:,:,count))
    colorbar
    drawnow

    count = count + 1;
end
% pulse_replica=m_chirp(Parameters.timeSeriesParameters(1)*.5,Parameters.timeSeriesParameters(2),0,Parameters.basebandSamplingRate);
% [h,i] = xcorr(pulse_replica);
% h = abs(h);
% [~,hc] = max(abs(h));
% h = h/norm(h);
% s = h(hc-length(I)/2:hc+length(I)/2-1);
% 
% psfb = kron(s,s');
% 
% for angle = [0 30 45 60 ]+5
%     p = imrotate(psfb,angle,'crop');
%     psf(:,:,count) =  p/sum(p(:));
%     B(:,:,count) = imfilter(I,psf(:,:,count)) ;
%     imagesc(B(:,:,count))
%     drawnow
% 
%     count = count + 1;
% end

g = B;
gn = B;
f = B(:,:,1);
fn = mean(B,3);
figure
imagesc(f)
 xlim([80 150])
    ylim([80 150])
    %%
    props.maxIter = 60;
    props.alpha = .7
    props.lambda = .04;
    props.P = 3;
    props.beta =.01;%

    HR=robustsr(B,fn,  psf, props);
    
%%
figure
for ii = 1:1000
    ii
    z = zeros(size(f));
    for jj = 1:size(psf,3)
        z = z + imfilter(g(:,:,jj)-gn(:,:,jj),psf(:,:,jj).^2);
    end
    fn = fn+  z;
    imagesc(fn)
    xlim([80 150])
    ylim([80 150])
    drawnow
    for jj = 1:size(psf,3)
        gn(:,:,jj) = imfilter(fn,psf(:,:,jj)) ;
    end

end