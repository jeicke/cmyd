rpos = [0 0 0];   
tpos = [1000 0 0];
c = [-500 0 0];
Fs = 60E6;
Fc = 9410E6;
% pulse rep frequency in HZ

rpm = 24;
bwh = 1.8;
prf = 3000;
Npulses = floor(1/(rpm/60) * bwh/360 * prf);
%dwell time in seconds
dwell_time = 1/(rpm/60) * bwh/360;

pulse_length = 1200E-9;
bw = 60E6;
radsim('C:\Users\rnl\Desktop\testdata\test',dwell_time,[0 0 0],[0 0 0],0,'propogate',@propagateDataFast,...
       'timeseries' ,@chirppulse,'transmitter_start' ,tpos,...
       'receiver_start',rpos,'parameters',[bw pulse_length prf],...
       'fs',Fs,'fc',Fc ,'clutter',c);  
%%
[data Fs Fc channels] = testread('C:\Users\rnl\Desktop\testdata\test');
%%
%FFTSIZE = 32768;
%[spectrum freqs time] = calculate_psd(data',Fs ,FFTSIZE,.5,hanning(FFTSIZE),false);
%imagesc(time,freqs/1E6,10*log10(abs(squeeze(spectrum))))

%%

 
pulse_replica=m_chirp(bw,pulse_length,0,Fs)'.';
pulse_replica=pulse_replica(end:-1:1);
 plength = round(Fs/ prf);
datapad = [data zeros(1,ceil(length(data)/ plength) *  plength-length(data))];
datapad = 100*datapad + 1/sqrt(2) * (randn(size(datapad)) + 1i * randn(size(datapad)));

D = reshape(datapad,plength,[]);
[D ] = pulsecompression(D,pulse_replica,'cheb',70);

[Df] = dopplerfilter(D,128,'cheb',70,false);
f = (prf/2)* (-size(Df,2)/2:size(Df,2)/2-1)/(size(Df,2)/2);
p = fftshift((20*log10(abs(Df))),2);
rngspa= 299792458/(Fs)/2;
ranges = [0:1:size(p,1)-1] * rngspa;
imagesc(f,ranges/1000,p)
caxis([max(p(:))-80 max(p(:))]);
colorbar
drawnow;