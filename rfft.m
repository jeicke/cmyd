function [x] = rfft(varargin)

if(length(varargin)>1)
    x = fft((varargin{1}),varargin{2});
else
    x = fft((varargin{1}));
end
x = x(1:length(x)/2+1,:);

