function [sinr sinr2] = calculate_sinr(ws,rsp,spatial_sv,Npulses,Nsensors,step,guard_bandwidth,snapshots,points_kept,gates,doppler_gate )
Wd = fftshift(dftmtx(Npulses),2);
M = Npulses-3+1;
u =   fftshift(dftmtx(M),2);
shade_doppler =  chebwin(M,70) ;
count = 1;
pulse = [0:Npulses-1]'/Npulses;

sinr = [];
pulse2 = [0:M-1]'/M;
L = ceil(snapshots/2);
rsp = permute(rsp,[1 3 2]);
for range_gate =gates 
    start_left = max([range_gate-(L+guard_bandwidth) 1 ]);
    end_left = max([range_gate - guard_bandwidth-1 1]);
    start_right = min([range_gate + guard_bandwidth+1 size(rsp,1)-1]);
    end_right = min([range_gate+(L+guard_bandwidth) size(rsp,1)-1]);
    s = (rsp([start_left:end_left],:,:));
    s2 = abs(s(:,:,floor(size(s,3)/2)));
    [mv iii]=  sort(max(s2,[],2));
    s = s(iii(end:-1:max([end-points_kept 1])),:,:);
    % s = s(iii(end),:,:);
    s = permute(reshape(s(:),size(s,1),Nsensors,Npulses,[]),[1 4 2 3]);
    s = reshape(s(:),[],Nsensors *Npulses);
    sl = s;
    s = (rsp([start_right:end_right],:,:));
    s2 = abs(s(:,:,floor(size(s,3)/2)));
    [mv iii]=  sort(max(s2,[],2));
    s = s(iii(end:-1:max([end-points_kept 1])),:,:);
    %s = s(iii(end),:,:);
    s = permute(reshape(s(:),size(s,1),Nsensors,Npulses,[]),[1 4 2 3]);
    s = reshape(s(:),[],Nsensors * Npulses);
    sr = s;
    s = [sl;sr];
    R = (s' * s)/size(s,1);
    
    W = zeros(Npulses * Nsensors,M);
    for ii = 1:M
        for jj = 1:3
            index_start = 1+(ii-1) * size(ws,1) + (jj-1) * size(ws,1);
            W(index_start:index_start + size(ws,1)-1,ii) = ws(:,jj,range_gate);
        end
    end
 %   e = squeeze(rsp(range_gate,:,:)).';
   % e = e(:);
   % R = e*e';
%    s = permute(reshape(rsp(:),size(rsp ,1),size(spatial_sv ,1),Npulses,[]),[1 4 2 3]);
% s = reshape(s(:),[],size(spatial_sv ,1) *Npulses);
% s = s(1000:1200,:);
% R = s'*s/size(s,1);
%R = R+1 * eye(size(R));
    R = R + eye(size(R));
% 
    Ri = inv(R);
    count2 = 1;

    %for doppler_gate = -200:200
    w = W *u;
    for ii = 1:size(Wd,2)
        %doppler_sv = exp(-1i * 2 * pi * doppler_gate/12 * pulse);
        doppler_sv = Wd(:,ii);
        sv = kron(doppler_sv,spatial_sv);%
        w2 = Ri * sv;

        %sinr( count) =  abs(w'*sv).^2/abs((w'*R*w));
        sinr( count,count2) =  abs(w2'*kron(spatial_sv,doppler_sv)).^2/abs(dot(w2,R*w2));
       %
%         for jj = 1:size(u,2)
%             uc = u(:,jj);
%             w =W * uc ;
% 
%            sinr_temp(jj) =  abs(w'*sv).^2/abs((w'*R*w));
%         end
        sinr_temp = abs(w'*sv).^2./abs(dot(w,R*w).');
        
        sinr2( count,count2) = max(sinr_temp);
        sinr( count) =  1*sinr( count)/(npulses * size(geometry,1));
        count2 = count2 + 1;
    end
%    for doppler_gate = -500:500
%         doppler_sv = exp(i * 2 * pi * doppler_gate/12 * pulse);
%         u = exp(i * 2 * pi * doppler_gate/12 * pulse2 );
%         sv = kron(spatial_sv,doppler_sv);%
% 
%     
% %         e = squeeze(interference(range_gate,:,:)).';
% %          e = e(:);
%          
%          w = W* u ;
%          w = Ri * sv;
%          
%          sinr(count,count2) = abs(w'*sv).^2/real(w'*R*w);%;%real((w'*s') * (s*w));%abs(w'*sv).^2;%/real((w'*s') * (s*w));%abs(w'*sv).^2;%
%        %  R = (W'*e)*(e'*W);
%         count2 = count2 +1;
%        % sinr(count,doppler_gate) =1;%(abs(f'*v)).^2/real(f'*R*f);%(abs(f'*v)).^2;%1/real(f'*R*f);%(abs(f'*v)).^2;%/;
%     end
    count = count + 1;
end