close all;
Fs = 20000000; %samples/second
Fc = 9350E6; %hertz
%                    b:    Linear FM Bandwidth (Hz)
%                    t:    Pulse width (seconds)
%                    f0:   Starting frequency (Hz)
%                    fs:   Sampling frequency (Hz)
h=m_chirp(5E6,6.4000e-06,-2.5E6,Fs);

x = zeros(1,4096);

x(1100) = 1;

%x(1107) = exp(1i * 2 * pi *rand(1));
%x(1100) = 0;


N = length(h)-1;
h = h/sqrt(h'*h);
for ii = 1:length(x)-N
    y(ii) = x(ii:ii+N)*h;
end
y = y + .00001*randn(size(y));
x = x + .00001;
%y = filter(h,1,s + .001*randn(size(s)));
m = filter(conj(h(end:-1:1)),1,y);

pc = m(900:1150);
X = x(900:1150);
plot(20*log10(abs(pc)));
fpc = fft(pc);
fpd = fft(X);
ei = 1i*2*pi/length(fpc)*[0:(length(fpc)-1)]';
e = exp(ei*([0:length(fpc)-1]))/length(fpc);
for ii= 1:size(e,2)
    e(:,ii) =  e(:,ii)/sqrt((e(:,ii)'*e(:,ii)));
end

pc2 = fpc*e;
K = zeros(length(fpc));
for ii = 1:100
    fpc = fpc + .0001*randn(size(fpc))  + 1i*.0001*randn(size(fpc));
    K = K + fpc'*fpc;
end
clear S S2;
opts.POSDEF = true; opts.SYM = true;
sigma = .02;
%  = (s' * s);
%         x = linsolve(R+ sigma * mean(diag(R))  * eye(size(R)),C,opts);
%         n =real(dot(C,x));
%         w =   x./ sqrt(repmat(n,size(x,1),1));
for ii= 1:size(e,2)
    e2 =  e(ii,:);
    
    [x] = linsolve(K+.001*eye(size(K)),e2);
    w = x./dot(e2,x);
    S(ii)= dot(w,K * w);

    %S(ii) = 1./(e2'*inv(K+.0001*eye(size(K)))*e2);
    
    S2(ii) = dot(e2,K*e2);
end
plot(10*log10(abs(S2)),'g');
hold on;
plot(10*log10(abs(S)),'r');