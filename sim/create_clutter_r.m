function [R, R2] = create_clutter_mimo(AZ,power,Dopplers,K,geometry,transmit_geometry,EL,...
    npulses,prf,range,unambiguous_range,power_threshold,transmitters,transmit_beam,range_ambiguous,range_gate,dopplerMax)

x = 10*log10(max(abs(power),[],2));
% trim off low power portions
[c ii] = find(x>(max(max(x))-power_threshold));

Tr = 1/prf;
pulse_sv=(0:(npulses-1))' * Tr;


%fc = transmitters/((npulses/transmitters) * 1/prf);

a0 = 1/(transmitters * Tr);

if(range_ambiguous)
    if(range_gate>length(unambiguous_range))
        range_gate = range_gate-length(unambiguous_range);

    end
   for ii = 1:ceil(size(power,2)/length(unambiguous_range))
       range_gate(ii+1) = range_gate(ii) + length(unambiguous_range);
   end
end
range_gate(find(range_gate>size(power,2))) = [];
%extra_doppler(1:1:end) = prf/2 * 2 * pi * 1i;
%prf = prf/transmitters;

Dopplers = (Dopplers(c,:));


power = (power(c,:));

voltage = (sqrt(abs(power)));%.*exp(i*2*pi*rand(size(power))));

A =sin(AZ(c));
B = cos(AZ(c));
C = cos(EL(range_gate));
D = sin(EL(range_gate));

S =  [kron(C,A);kron(C,B);kron(ones(size(A)),D)];
SV = exp(double(-2 * pi * K * geometry * S * 1i));
B = reshape(SV,[],length(c),length(EL(range_gate)));
if(~isempty(transmit_geometry))
    SV = exp(-2 * pi * K * squeeze(transmit_geometry)' * S * 1i);

    C = reshape(SV,[],length(c),length(EL(range_gate)));
end
transmit_beam = transmit_beam * pi/180;
% extra_doppler = zeros(size(C,2),size(C,3),size(pulse_sv));
% for ii = 1:transmitters
%     l = mod(ii,transmitters);
%     SV_transmit = squeeze(C(ii,:,:));
%     
% 
%     extra_doppler = extra_doppler + SV_transmit .* exp(j * 2 * pi * pulse_sv * l * a0);
% end



Dopplers = Dopplers(:,range_gate);
if(~isempty(dopplerMax))
dopplerMax = dopplerMax(range_gate);
end
voltage = double(voltage(:,range_gate));
voltage = voltage;
clear R;

for jj = 1:length(range_gate)

    for ii = 1:npulses
        extra_doppler = zeros(size(Dopplers(:,jj)));
        if(transmitters>1)
            for kk = 1:transmitters
                l = mod(kk,transmitters);
                SV_transmit = squeeze(C(kk,:,jj)).';
                extra_doppler = extra_doppler + SV_transmit .* exp(j * 2 * pi * pulse_sv(ii) * l * a0);
            end
        else
            extra_doppler= 1;
        end
        extra_doppler =  extra_doppler / transmitters;
        if(isempty(dopplerMax)||numel(range_gate)==1)
            doppler(ii,:) = voltage(:,jj).*extra_doppler.*exp(2 * pi * 1i * pulse_sv(ii)*(Dopplers(:,jj)));
        else
            refGate = round(length(range_gate)/2);
            crv  = Dopplers(:,jj)-(dopplerMax(jj)-dopplerMax(refGate));
            doppler(ii,:) = voltage(:,jj).*extra_doppler.*exp(2 * pi * 1i * pulse_sv(ii)*(crv));
        end
        
        %((SV*( voltage.* exp(pulse_sv(ii)*Dopplers)))).';%((SV*( voltage.* exp(pulse_sv(ii)*Dopplers)))).';%((SV*( voltage.* exp(pulse_sv(ii)*Dopplers))).*SV_el).';
    end
    clear s;
    for ii = 1:size(doppler,2)
      s(:,ii)  =  kron(B(:,ii,jj),doppler(:,ii));
    end
    R(:,:,jj) = s ;
    for ii = 1:size(doppler,2)
      s(:,ii)  =  kron(doppler(:,ii),B(:,ii,jj));
    end
    R2(:,:,jj) = s;
end
% R2 = R;
% clear R;
% for jj =1:length(range_gate)
%     R(:,:,jj) = R2(:,:,jj) * R2(:,:,jj)';
% end