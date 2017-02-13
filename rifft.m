function [x] = rifft(y)
y = [y; y(end-1:-1:2,:)'.'];

x = ifft(y);