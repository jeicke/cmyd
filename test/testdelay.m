
Fs = 1024;
L = 64;

s = zeros(L,1);
s(1:7) = hann(7);
ss = fft(s);
s2 = rfft(s);
wave_number = single(-1i * 2 * pi  *[0:1:L/2-1 0] * (Fs/L)).';
for delay = 0:.01:16
     H1 = exp( wave_number* delay/Fs);
     h = sinc([0:L-1]-delay)' ;
     H = fft(h);
     m = ss.*H;
     mnf = s2.*H1;
     
     mmnf = rifft(mnf);
     mm = ifft(m);
     hold off;
     
     plot(real(mm),'bx-')
     hold on;
     plot(real(mmnf),'co-')
     plot(1000*imag(mm),'gx-')
     hold on;
     plot(s,'r')
     axis([1 16 0 1.5])
     
     pause(.01)
     if(delay == floor(delay))
         pause(1)
     end
end