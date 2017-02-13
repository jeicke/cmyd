function [sinr] = sinrR_postdoppler(R,Npulses,Nsensors,gates,V,prf)
pulse=((0:(Npulses-1))')/prf;
z = repmat(pulse,1,length(gates));
g = repmat(gates,length(pulse),1);
doppler_sv = exp(-1i * 2 * pi * z.*g);
window = repmat(chebwin(size(doppler_sv ,1),80),1,size(doppler_sv ,2)) ;
doppler_sv  = doppler_sv .* window;
R2 = reshape(R(:),Npulses,[]);
pdR = doppler_sv.' * R2;
opts.POSDEF = true; opts.SYM = true;
for ii = 1:size(pdR,1)
    s = pdR(ii,:);
    s = reshape(s,Nsensors,[]);
    R = s * s'/size(s,2);
    noise = .001;
    R = R+ noise * eye(size(R));
    w = linsolve(R,V,opts);

    sinr(ii)=  noise  * abs(dot(w,V)).^2./abs(dot(w,R*w));
    
end
sinr = sinr/( Nsensors);
