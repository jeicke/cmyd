function [sinr,w] = sinrR(R,Npulses,Nsensors,gates,V,prf)
pulse=((0:(Npulses-1))')/prf;
R = sum((R),3);
R = double(R) ;

sigma =.001;
noise =sigma * mean(diag(R));
R = R+noise* eye(size(R)) ;
%R = Rn;
%Ri = inv(R);
count = 1;
opts.POSDEF = true; opts.SYM = true;
if(~isempty(gates))
    z = repmat(pulse,1,length(gates));
    g = repmat(gates,length(pulse),1);
    doppler_sv = exp(1i * 2 * pi * z.*g);
%     count = 1;
%     for dv = -prf:prf
%         doppler_sv(:,count) = exp(1i*2*pi*dv*([1:Npulses]-1)/prf);
%         count = count + 1;
%     end
    sv = kron(V,doppler_sv);%
   
    w = linsolve(R,sv,opts);
    %w = w./(dot(w,R*w));
    sinr=  noise*abs(dot(w,sv)).^2./abs(dot(w,R*w));
    sinr =  sinr/(Npulses * Nsensors);
   
        
%     clear sinr;
%     for doppler_gate = gates
%         doppler_sv = exp(i * 2 * pi * doppler_gate * pulse);
%         if(~exist('V','var')||isempty(V))
%             sv = kron(ones(Nsensors,1) ,doppler_sv);%
%         else
%             sv = kron(V ,doppler_sv);%
%         end
%         w = Ri * sv;
% 
%         %sinr( count) =  abs(w'*sv).^2/abs((w'*R*w));
%         sinr( count) =  noise  * abs(w'*sv).^2/abs((w'*R*w));
%         sinr( count) =  sinr( count)/(Npulses * Nsensors);
%         count = count + 1;
%     end
% a = 1;
else
   
    R = double(R) + noise * eye(size(R));
    Ri = inv(R);
    count = 1;
    Wd= fftshift(dftmtx(Npulses),2);
    for ii = 1:size(Wd,2)
        %doppler_sv = exp(i * 2 * pi * doppler_gate/13.333* pulse);
  
        if(~exist('V','var')||isempty(V))
            sv = kron(ones(Nsensors,1) ,Wd(:,ii));%
        else
            sv = kron(V ,Wd(:,ii));%
        end
        w = Ri * sv;

        %sinr( count) =  abs(w'*sv).^2/abs((w'*R*w));
        sinr( count) =  abs(w'*sv).^2/abs((w'*R*w));
        sinr( count) = noise*sinr( count)/(Npulses * Nsensors);
        wx (:,ii)= w;
        count = count + 1;
    end
    w = wx;
    
end
