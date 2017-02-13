function [sinr,w] = sinrRpdm(R,Npulses,Nsensors,gates,V,prf,pulses)
x =Npulses-pulses+1;
pulse=((0:(x-1))')/prf;

x =Npulses;
pulse2=((0:(x-1))')/prf;
R = sum((R),3);
R = double(R) ;

sigma =.002;
noise =sigma * mean(diag(R));
R = R+noise* eye(size(R)) ;
%R = Rn;
%Ri = inv(R);
count = 1;
opts.POSDEF = true; opts.SYM = true;

z = repmat(pulse,1,length(gates));
g = repmat(gates,length(pulse),1);
doppler_sv = exp(1i * 2 * pi * z.*g);
window = repmat(chebwin(size(doppler_sv ,1),70),1,size(doppler_sv ,2)) ;
doppler_sv  = doppler_sv .* window;

z = repmat(pulse2,1,length(gates));
g = repmat(gates,length(pulse2),1);
doppler_sv2 = exp(1i * 2 * pi * z.*g);
window = repmat(chebwin(size(doppler_sv2 ,1),70),1,size(doppler_sv ,2)) ;
doppler_sv2  = doppler_sv2 .* window;

In = eye(Nsensors);
for d = 1:size(doppler_sv,2)
    Fm = zeros(size(doppler_sv,1),pulses);
    for pdb =  1:pulses
        
        endIndex = Npulses-(pulses-pdb);
        Fm(pdb:endIndex,pdb) = doppler_sv(:,d);
    end
    
    T = kron(In/6,Fm);
    Rd = T'*R*T;
    sv = kron(V,doppler_sv2(:,d) );%
    w = linsolve(Rd,T'*sv,opts);
    
    w = T*w;
    sinr(d)=  abs(dot(w,sv)).^2./abs(dot(w,R*w));
    sinr(d) =  sinr(d)/(Npulses * Nsensors);
end
a = 1;




