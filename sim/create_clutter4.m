% power is range X azimuth matrix of 
function [clutter] = create_clutter4(AZ,power,Dopplers,K,geometry,EL,npulses,prf,range,unambiguous_range,power_threshold,range_ambiguous)

x = 10*log10(max(abs(power),[],2));
% trim off low power portions
[c ii] = find(x>(max(max(x))-power_threshold));

pulse_sv=single((0:(npulses-1))'/prf);

Dopplers = (Dopplers(c,:));

Dopplers = single(Dopplers * 2 * pi) * 1i;
power = (power(c,:));





voltage = (sqrt(abs(power)).*exp(i*2*pi*rand(size(power))));

A =sin(AZ(c));
B = cos(AZ(c));
C = cos(EL);
D = sin(EL);

SV = exp(- 2 * pi * K * geometry * [kron(C,A);kron(C,B);kron(ones(size(A)),D)]* 1i);
B = reshape(SV,[],length(c),length(EL));
clutter = single(zeros(size(Dopplers,2),npulses,size(geometry,1)));

for jj = 1:size(Dopplers,2)

    for ii = 1:npulses
        clutter(jj,ii,:) =  B(:,:,jj) * (voltage(:,jj).* exp(pulse_sv(ii)*Dopplers(:,jj)));
        %((SV*( voltage.* exp(pulse_sv(ii)*Dopplers)))).';%((SV*( voltage.* exp(pulse_sv(ii)*Dopplers)))).';%((SV*( voltage.* exp(pulse_sv(ii)*Dopplers))).*SV_el).';
    end
end
if(range_ambiguous)
clutter_reduced = zeros(length(unambiguous_range),size(clutter,2),size(clutter,3));
for ii = 1:ceil(size(clutter,1)/length(unambiguous_range))
    index = (1:length(unambiguous_range)) + (ii-1) * length(unambiguous_range);
    %index = mod(index,length(unambiguous_range)) + 1;
    index(find(index>size(clutter,1)))= [];
    index2 = 1:(min([length(range) length(index)]));
    clutter_reduced(index2,:,:) = clutter_reduced(index2,:,:)  + clutter(index,:,:);
end
clutter = clutter_reduced;
end
% a = 1;
%now find maximum unambiguous range index
% s = permute(reshape(clutter(:),size(clutter ,1),size(SV ,1),npulses,[]),[1 4 2 3]);
% s = reshape(s(:),[],size(SV ,1) *npulses);
% s = s(1000:1200,:);
% R = s'*s/size(s,1);
% R = R+1 * eye(size(R));
% 
% Ri = inv(R);
% 
% 
% pulse = [0:npulses-1]'/npulses;
% % SVa =2 * pi * K * [sin(AZ(39));cos(AZ(39));0];
% % 
% % SVa = exp(geometry * SVa * 1i);
% % SVa = SV(:,39)'.';
% % count = 1;
% % doppler_sv = exp(i * 2 * pi * 0/12 * pulse);
%  
% % for az = [-90:.25:90]
% %     SVa =2 * pi * K * [sin(az * pi/180)*cos(el * pi/180);cos(az*pi/180)*cos(el * pi/180);sin(el*pi/180)];
% % 
% %     SVa = exp(geometry * SVa * 1i);
% %     
% %     sv = kron(SVa ,doppler_sv);%
% %     w = Ri * sv;
% % 
% %     %sinr( count) =  abs(w'*sv).^2/abs((w'*R*w));
% %     sinr( count) =  sv' * R * sv;%abs(w'*sv);%sv' * Ri * sv;%abs(w'*sv).^2/abs((w'*R*w));
% %     sinr( count) =  1*sinr( count)/(npulses * size(geometry,1));
% %     count = count + 1;
% % end
% el = EL(1100);
% az = 0;
% SVa =2 * pi * K * [sin(az * pi/180)*cos(el );cos(az*pi/180)*cos(el );sin(el)];
% SVa = exp(geometry * SVa * 1i);
% 
% count = 1;
% for doppler_gate = -200:200
%     doppler_sv = exp(i * 2 * pi * doppler_gate/12 * pulse);
%     sv = kron(SVa ,doppler_sv);%
%     w = Ri * sv;
% 
%     %sinr( count) =  abs(w'*sv).^2/abs((w'*R*w));
%     sinr( count) =  abs(w'*sv).^2/abs((w'*R*w));
%     %sinr( count) =  1*sinr( count)/(npulses * size(geometry,1));
%     count = count + 1;
% end
% a = 1;


% clutter = zeros(size(power,2),npulses,size(SV ,1));
% for ii = 1:size(SV ,1)
%     temp = repmat(squeeze(SV(ii,:)),[],size(power,2)).';
%     temp_ch =  sum((power .* temp),1);
% 
%     clutter(:,:,ii) = clutter_pulse.*repmat(temp_ch,npulses,[]).' ;
% end
% a = 1;