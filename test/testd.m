

a = [1 2 3 2 1];


b = zeros(128,1);
b(1) = 1;
aa = zeros(size(b))

aa(1:length(a)) = a;
aa = circshift(aa,-2);
c = ifft(fft(aa).*fft(b));
plot(c);