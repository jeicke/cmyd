function resp = calcBP(posit,wts,steer)
theta = [-90:0.1:90]*pi/180;
posit = posit*2*pi;
e = exp(i*posit*sin(theta));
wts = wts.*exp(i*posit*sin(steer));
resp = wts'*e;
resp = resp.*sqrt(cos(theta));
%resp = abs(resp).^2;
%figure
%plot(theta,20*log10(abs(resp)))
%axis([-1.75 1.75 -80 5]);
