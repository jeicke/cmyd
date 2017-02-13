
function clutter = create_clutter3(AZ,power,Dopplers,K,geometry,EL,npulses,prf,range,unambiguous_range,power_threshold)


A =sin(AZ);
B = cos(AZ);
C = cos(EL);
D = sin(EL);
SV =2 * pi * K * [A;B;zeros(size(A))];

SV = exp(geometry * SV * 1i);

SV_el = single(exp(geometry * [C;C;D] * 2 * pi * K * 1i)); 
pulse_sv=single((0:(npulses-1))'/prf);

clutter_pulse = zeros(size(power,2),npulses);






Dopplers = single(Dopplers * 2 * pi) * 1i;
x = 10*log10(max(abs(power),[],2));
% trim off low power portions
[c ii] = find(x>(max(max(x))-power_threshold));
power = power(c,:);

Dopplers = Dopplers(c,:);

SV = SV(:,c);

voltage = sqrt((power)) .* exp(i* 2*pi * rand(size(power)));
clutter = single(zeros(length(c),npulses,size(SV ,1)));
for ii = 1:npulses
     e =voltage.* exp(pulse_sv(ii)*Dopplers);
     e = repmat(e.',size(SV,1),[]);
     clutter(:,ii,:) = (SV.*e).';%voltage.* exp(pulse_sv(ii)*Dopplers);
end
s = permute(reshape(clutter(:),size(clutter ,1),size(SV ,1),npulses,[]),[1 4 2 3]);
s = reshape(s(:),[],size(SV ,1) *npulses);

R = s'*s/size(s,1);
R = R+1 * eye(size(R));

Ri = inv(R);

AZ(51) = pi/2;
pulse = [0:npulses-1]'/npulses;
EL(500) = 0;
SVa =2 * pi * K * [1;0;0];

SVa = exp(geometry * SVa * 1i);
SVa = SV(:,78);
count = 1;
D = fftshift(dftmtx(npulses),1);
for doppler_gate = 1:size(D,1)
    doppler_sv = D(doppler_gate,:).';%exp(i * 2 * pi * doppler_gate/12 * pulse);
    sv = kron(SVa,doppler_sv);%
    w = Ri * sv;

    sinr( count) =  sv' * Ri * sv;%abs(w'*sv).^2/abs((w'*R*w));
    sinr( count) =  1*sinr( count)/(npulses * size(geometry,1));
    count = count + 1;
end
a = 1;
count = 1;
for doppler_gate = -200:200
    doppler_sv = exp(i * 2 * pi * doppler_gate/12 * pulse);
    sv = kron(SVa ,doppler_sv);%
    %w = Ri * sv;

    %sinr( count) =  abs(w'*sv).^2/abs((w'*R*w));
    sinr( count) =  sv' * Ri * sv;%abs(w'*sv).^2/abs((w'*R*w));
    sinr( count) =  1*sinr( count)/(npulses * size(geometry,1));
    count = count + 1;
end
a = 1;
% for ii = 1:npulses
% 
%     clutter(:,ii,:) = ((SV*( voltage.* exp(pulse_sv(ii)*Dopplers))).*SV_el).';
% 
% end
%now find maximum unambiguous range index


