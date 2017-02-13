function [wgt] = chebwgt(nel,slob)
% Function chebwgt calculates the Chebyshev weights using the chirp-z
% transform in Function czt1ds. The sidelobe level is input in dB in
% slob. The weights can be of either odd or even length.
% Inputs:
%	 nel : desired length
%	slob : sidelobe level
% Output:
%	 wgt = Chebyshev weights normalized to unit energy
% This function call czt1ds
% Format wgt = chebwgt(nel,slob)

xn1 = nel-1 ;
nbar = .5*xn1 ;
nfac = 1./nel ;
wr = 10.^(abs(slob)/20.) ;
ws = wr+sqrt(wr*wr-1.0) ;
ws=log(ws)/xn1 ;
ws = exp(ws) ;
mu = 0.5 * (ws+1.0/ws) ;
cp0 = 1.0/wr ;
qa = pi*nfac ;
w = mu*cos((0:xn1)*qa) ;

 for i = 1:nel
    if(abs(w(i))<=1.0)
       array(i) = cp0 * cos(xn1*acos(w(i))) ;
    else
       array(i) = cp0 * cosh(xn1*acosh(w(i))) ;
    end
 end
 u2 = nbar*nfac ;
 u1 = -u2 ;
 array = czt1ds(array,nel,nel,u1,u2,1,1) ;
 wgt = real(array)' ;
 temp=sum(wgt)^2;
 wgt=wgt/sqrt(temp);
 %norm = 1/sqrt(sum(wgt.^2)) ;
 %wgt = norm * wgt ;
 
% plot weight response
%farray = array ;
%farray(nel+1:256) = zeros(1,256-nel) ;
%fwgt = 10*log10(abs(fft(farray)).^2+1.e-7);
%fwgt = fftshift(fwgt) ;
%plot((0:255),fwgt);

