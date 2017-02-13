% ==============================================================================
%
%                           ===>>>  CZT1DS.M  <<<===
%_______________________________________________________________________________
%
%   PROGRAMMER : Jim Ward
%
%   DATE   CODE: 29 August 91
%   UPDATE CODE: 29 August 91
%_______________________________________________________________________________
%
%   DESCRIPTION: Function czt1ds calculates Chebyshev weights for either odd or
%                even number using the chirp-z transform. This program has been
%                adapted from a program by E. Kelly running in the mainframe
%                at the Lincoln Laboratory. The program is used in conjunction
%                with function chebwgt. Function czt1ds simply does the chirp-z
%                transform.
%
%_______________________________________________________________________________
%
%   USAGE: wgt = czt1ds(array,nx,nupts,umin,umax,beta,isign)
%_______________________________________________________________________________
%
%   INPUTS : array    : Input sequence.
%            nx       : Input length.
%            nupts    : Output length.
%            umin     : Beginning of chirp-z transform.
%            umax     : Ending of chirp-z transform.
%            beta     : 1
%            isign    : 1
%
%   OUTPUTS: wgt      : Output sequence of length nupts.
%_______________________________________________________________________________
%
%   CALLED BY   : chebwgt
%
%   CALLS TO    : NONE
%
% ==============================================================================

function wgt = czt1ds(array,nx,nupts,umin,umax,beta,isign)

j = sqrt(-1) ;
mx = max(nx,nupts) ;

if (nupts==1)
   uinc = 0 ;
else
   uinc = (umax-umin)/(nupts-1) ;
end

lupts = 1 ;
logu  = 0 ;

while (lupts<(nx+nupts-1))
    lupts = 2 * lupts ;
    logu = logu + 1 ;
end

zbuf = zeros(size(1:lupts)) ;

w0 = isign * pi * uinc * beta ;
wu = w0 * (0:mx-1).^2 ;

w0 = 1 / lupts ;

arg(1:nupts) = -wu(1:nupts) ;
wt(1:nupts)  = w0 * ones(1,nupts) ;

arg(nupts+1:lupts-nx+1) = zbuf(nupts+1:lupts-nx+1) ;
wt(nupts+1:lupts-nx+1)  = zbuf(nupts+1:lupts-nx+1) ;

arg(lupts-nx+2:lupts) = -wu(nx:-1:2) ;
wt(lupts-nx+2:lupts)  = w0 * ones(1,nx-1) ;

wuray = wt .* exp(j*arg) ;
wuray = fft(wuray) ;

w0   = isign * 2. * pi * umin * beta ;
data = array .* exp(j*(w0*(0:nx-1)+wu)) ;
data(nx+1:lupts) = zbuf(nx+1:lupts) ;
data = fft(data) ;
data = data .* wuray ;

data = ifft(data) ;
wgt  = data(1:nupts) .* exp(j*wu(1:nupts)) ;

return
