function [sinr,dopplerbins,R] = theorysinr(track,geometry,prf,npulses,lamda,clutterarea,tvelocity,rvelocity)
[r_transmitter_scatter, aztransmit, detransmit]= computeRange(track, 0,1);
r_transmitter_scatter = r_transmitter_scatter(2:end);
% next comput range from scatterer to receiver and the angles
%bsxfun(@plus,scatterPos,-antPos)
[r_scatter_ant, azreceive, dereceive]= computeRange(track, 0,2);
r_scatter_ant = r_scatter_ant(2:end);
voltage =clutterarea.'./(r_transmitter_scatter.*r_scatter_ant);

azreceive = azreceive(2:end);
dereceive = dereceive(2:end);
aztransmit = aztransmit(2:end);
detransmit = detransmit(2:end);
A = sin(azreceive);
B = cos(azreceive);
C = cos(dereceive);
D = sin(dereceive);

At = sin(aztransmit);
Bt = cos(aztransmit);
Ct = cos(detransmit);
Dt = sin(detransmit);


clutterDoppler=   [A.*C B.*C  D] *rvelocity'* 1/lamda+[At.*Ct Bt.*Ct  Dt] *tvelocity'* 1/lamda;

S =  [A.*C B.*C D];
SV = exp(double(-2 * pi * 1/lamda *  S * geometry'*1i));
Tr = 1/prf;
pulseSv=(0:(npulses-1))' * Tr;
doppler = exp(2 * pi * 1i * pulseSv*clutterDoppler.');
clear s;
for ii = 1:size(doppler,2)
    s(:,ii)  =  kron(SV(ii,:).',doppler(:,ii));
end
R = s*s';
R = R/size(s,2);
dopplerbins = [-prf/2:1:prf/2];
%dopplerbins = linspace(-prf/2,prf/2,npulses);
[sinr{1}] = sinrR(R,npulses,size(SV,2),dopplerbins,ones(size(SV,2),1),prf);
[sinr{2}] = sinrRc(R,npulses,size(SV,2),dopplerbins,ones(size(SV,2),1),prf);
[sinr{3}] = sinrRpdm2(R,npulses,size(SV,2),dopplerbins,ones(size(SV,2),1),prf,1);
[sinr{4}] = sinrRpdm2(R,npulses,size(SV,2),dopplerbins,ones(size(SV,2),1),prf,2);
for ii = 1:size(doppler,2)
    s(:,ii)  =  kron(doppler(:,ii),SV(ii,:).');
end
R = s*s';
R = R/size(s,2);
