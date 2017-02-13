Fs = 20000000; %samples/second
Fc = 9350E6; %hertz
%                    b:    Linear FM Bandwidth (Hz)
%                    t:    Pulse width (seconds)
%                    f0:   Starting frequency (Hz)
%                    fs:   Sampling frequency (Hz)
h=m_chirp(5E6,6.4000e-06,-2.5E6,Fs);

x = zeros(1,4096);

x(1024) = 1;
x(1030) = 1;
x(1100) = 0;


N = length(h)-1;
for ii = 1:length(x)-N
    y(ii) = x(ii:ii+N)*h;
end
y = y + .001*randn(size(y));
x = x + .001;
%y = filter(h,1,s + .001*randn(size(s)));
m = filter(conj(h(end:-1:1)),1,y);

a = [zeros(1,N-2) h.' zeros(1,N-2)];
count = 1;
e = zeros(N+1,253);
for n = -N+1:N-1
    z = circshift(a.',n).';
    e(:,count) = z(N:2*N);
    count = count + 1;
end
count = 1;
clear m;
index = [1:(length(y)-N)];
 w =h/(h'*h);
for ii = index
   
    m(count) = w' * y(ii:ii+N).';
    
    count = count + 1;
end

clear mm;
clear mmm;
count = 1;
for ii = 700:1200
    
    p =x([-N+1:N-1]+ii+(length(h)-1)).';
    %p =m([-N+1:N-1]+ii).';
    rho = diag(abs(p).^2);
    
    R = e*rho*e';
    R = R + .0001 * mean(diag(R)) + eye(size(R));
    %R = diag(ones(size(h)));
    z = R\h;
    w =z./dot(h,z);
    mm(count) = w' * y(ii:ii+N).';
    
    mmm(count) = m(ii);
    count = count + 1;
end


%close all;
plot(20*log10(abs(mm)),'c')
hold on;
plot(20*log10(abs(mmm)),'r')
plot(20*log10(eps+abs(x([700:1200]+length(h)))),'b')

