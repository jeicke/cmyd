function [sinr,w] = sinrRc(R,Npulses,Nsensors,gates,V,prf)
pulse=((0:(Npulses-1))')/prf;
R = sum((R),3);
R = double(R) ;

sigma =.02;
noise =sigma * mean(diag(R));
R = R+noise* eye(size(R)) ;
%R = Rn;
%Ri = inv(R);
count = 1;
opts.POSDEF = true; opts.SYM = true;

z = repmat(pulse,1,length(gates));
g = repmat(gates,length(pulse),1);
windows =1;% chebwin(Npulses,80);
doppler_sv = exp(1i * 2 * pi * z.*g);
%     count = 1;
%     for dv = -prf:prf
%         doppler_sv(:,count) = exp(1i*2*pi*dv*([1:Npulses]-1)/prf);
%         count = count + 1;
%     end
sv = kron(V,bsxfun(@times,doppler_sv,windows));%

w = linsolve(eye(size(R)),sv,opts);
w = w./sqrt(dot(w,R*w));
sinr=  noise*abs(dot(w,sv)).^2./abs(dot(w,R*w));
sinr =  sinr/(Npulses * Nsensors);



