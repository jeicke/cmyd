% power is range X azimuth matrix of 
function [power Dopplers SV SV_el] = create_clutter2(AZ,power,Dopplers,K,geometry,EL,npulses,prf,range,unambiguous_range,power_threshold)
% trim off low power portions
x = 10*log10(max(abs(power),[],2));
[c ii] = find(x>(max(max(x))-power_threshold));
power = power(c,:);
Dopplers = Dopplers(c,:);
AZ = AZ(c);

SV =2 * pi * K * [sin(AZ);cos(AZ);zeros(size(AZ))];

SV = exp(geometry * SV * 1i);

SV_el = single(exp(geometry * [cos(EL);cos(EL);sin(EL)] * 2 * pi * K * 1i)); 


clutter_pulse = zeros(size(power,2),npulses);






Dopplers = single(Dopplers * 2 * pi) * 1i;


[c ii] = find(x>(max(max(x))-power_threshold));
power = power(c,:);
Dopplers = Dopplers(c,:);

SV = SV(:,c);
voltage = sqrt((power));
clutter = single(zeros(length(c),npulses,size(SV ,1)));
for ii = 1:npulses
     e =voltage.* exp(pulse_sv(ii)*Dopplers);
     e = repmat(e.',size(SV,1),[]);
     clutter(:,ii,:) = (SV.*e).';%voltage.* exp(pulse_sv(ii)*Dopplers);
end
% s = permute(reshape(clutter(:),size(clutter ,1),size(SV ,1),npulses,[]),[1 4 2 3]);
% s = reshape(s(:),[],size(SV ,1) *npulses);
% 
% R = s'*s/size(s,1);
% R = R+eye(size(R));
% 
% Ri = inv(R);
% 
% 
% pulse = [0:npulses-1]'/npulses;
% SVa =2 * pi * K * [0;1;0];
% 
% SVa = exp(geometry * SVa * 1i);
% count = 1;
% for doppler_gate = -100:100
%     doppler_sv = exp(i * 2 * pi * doppler_gate/12 * pulse);
%     sv = kron(SVa ,doppler_sv);%
%     w = Ri * sv;
% 
%     sinr( count) =  abs(w'*sv).^2/abs((w'*R*w));
%     count = count + 1;
% end
% a = 1;
% for ii = 1:npulses
% 
%     clutter(:,ii,:) = ((SV*( voltage.* exp(pulse_sv(ii)*Dopplers))).*SV_el).';
% 
% end
%now find maximum unambiguous range index
% power is range X azimuth matrix of 