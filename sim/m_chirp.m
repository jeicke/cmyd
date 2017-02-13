function y=chirp(b,t,f0,fs)
%
%       CHIRP - Compute Linear FM Pulse Data
%
%               Calling sequence:
%
%                    y = chirp( b , t , f0 , fs )
%
%               Inputs:
%
%                    b:    Linear FM Bandwidth (Hz)
%                    t:    Pulse width (seconds)
%                    f0:   Starting frequency (Hz)
%                    fs:   Sampling frequency (Hz)
%
%               Outputs:
%
%                    y:    Linear FM Pulse (complex sequence)
%
npc=fs*t;
w0=2*pi*f0/fs;
dw=2*pi*b/fs;
mu=dw/npc;
v=(0:npc-1)';
phi=w0*v + (mu/2)*v.^2;
y=exp(j*phi);
