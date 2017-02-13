function [y] = MakeChirp(bw, len, Fs, dir)
if(dir ==0)
f0 = -(bw/2);
fband = bw;

dur = len/Fs;

t = (0:(len-1))'/Fs;
phi = f0*t + (fband/(2*dur))*t.^2;
y = exp(2*pi*1i*phi);

% normalize the pulse %
y = y / norm(y);
else
    f0 = (bw/2);
fband = bw;

dur = len/Fs;

t = (0:(len-1))'/Fs;
phi = f0*t - (fband/(2*dur))*t.^2;
y = exp(2*pi*1i*phi);

% normalize the pulse %
y = y / norm(y);
end