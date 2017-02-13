function [sinr] = sinr_predoppler(rsp,Npulses,Nsensors,gates,V,Vdoppler,prf,Rall)
C = zeros(size(rsp,3) * size(Vdoppler,1),size(Vdoppler,2) * size(V,2));
Rall = sum(Rall,3);
count = 1;
for jj = 1:size(Vdoppler,2)
    for ii = 1:size(V,2)
        C(:,count) = colkron(Vdoppler(:,jj),V(:,ii));
        count = count + 1;
    end
end

opts.POSDEF = true; opts.SYM = true;

rsp2 = rsp;
R = zeros(size(Vdoppler,1) *Nsensors,size(Vdoppler,1) *Nsensors);
for ii = 1:size(rsp2,4)
    rsp = permute(rsp2(:,:,:,ii),[1 3 2]);


    s = (rsp(:,:,:));
    s = permute(reshape(s(:),size(s,1),size(V,1),size(Vdoppler,1),[]),[1 4 2 3]);
    s = reshape(s(:),[],size(V,1) * size(Vdoppler,1));
    R = R + (s' * s)/size(s,1);
end
noise  = 	100;
x = linsolve(R+ noise * eye(size(R)),C,opts);
n =real(dot(C,x));
ws=   x./ sqrt(repmat(n,size(x,1),1));

ws = reshape(ws(:),Nsensors,size(ws,1)/Nsensors);

M = Npulses-size(Vdoppler,1)+1;
W = zeros(Npulses * Nsensors,M);

for jj = 1:size(Vdoppler,1);
for ii = 1:M

        index_start = 1+(ii-1) * Nsensors + (jj-1) * Nsensors;
        W(index_start:index_start + Nsensors-1,ii) = ws(:,jj );
        

    end
end

u =  dftmtx(M);

pulse=((0:(Npulses-1))')/prf;
pulse2=((0:(Npulses-size(Vdoppler,1)))')/prf;
count = 1;

%w = W * u;
Wd = fftshift(dftmtx(Npulses),2);

z = repmat(pulse,1,length(gates));
g = repmat(gates,length(pulse),1);
doppler_sv = exp(-i * 2 * pi * z.*g);
z = repmat(pulse2,1,length(gates));
g = repmat(gates,length(pulse2),1);
doppler_sv2 = exp(-i * 2 * pi * z.*g);
w = W * doppler_sv2 ;
sv = kron(doppler_sv,V);%
sinr_temp = noise * abs(dot(w,sv)).^2./ abs(dot(w'.',(Rall +noise* eye(size(Rall)))*w'.'));
sinr = sinr_temp/(Npulses * Nsensors);


% w = W * doppler_sv2 ;
% tic
% for doppler_gate = gates
% 
%     doppler_sv = exp(-i * 2 * pi * doppler_gate * pulse);% .*chebwin(length(pulse),40) ;
%     doppler_sv2 = exp(-i * 2 * pi * doppler_gate * pulse2);
% 
%     w = W * doppler_sv2 ;
%     sv = kron(doppler_sv,V);%
% 
%    % sinr_temp = abs(sv' * Rall * sv);
%     
%    sinr_temp = abs(w'*sv).^2./ abs(dot(w'.',(Rall + eye(size(Rall)))*w'.')).';%./abs(w'*(Rall + eye(size(Rall)))*w);
%    %sinr_temp = abs(dot(w,(Rall + eye(size(Rall)))*w)).';
%         
%     sinr(count) = max(sinr_temp)/(Npulses * Nsensors);
%     count = count + 1;
%    
% end
% toc
% a = 1;
