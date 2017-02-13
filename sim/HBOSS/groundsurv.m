function [eangles,gangles,sranges]=groundsurv(height,minr,maxr,bw);
%distances in km, elevation beamwidth in sin(theta)

iter=1;
sranges(iter)=maxr;
[e,g]=sargmtiangles(height,sranges(iter));
eangles(iter)=e; 
gangles(iter)=g;
while sranges(iter)>minr
    iter=iter+1;
    eangles(iter)=asin(sin(e*pi/180)-bw)*180/pi;
    sranges(iter)=sargmtiangles2(height,eangles(iter));
    [e,g]=sargmtiangles(height,sranges(iter));
    gangles(iter)=g;
end

