function target = target_generate_mimo(target_position,platform_position,velocity,...
                                  pulse_replica,ranges,BW,beam_center,...
                                  kT,losses,T,NF,rcs,...
                                  Ptransmit,PSG_transmit,frequency,npulses,...
                                  prf,transmit_geometry,fullarray_geometry,...
                                  subarray_geometry,subarray_window,mimo_transmit_geometry,psg_subarray,unambiguous_range,transmitters)
wavelength = 299792458/frequency;
rngspa= 299792458/(BW)/2;


range = 1000*sqrt(sum((target_position-platform_position).^2,1));
geometry= fullarray_geometry;
clear HGA;
HGA.theta_source = beam_center(1);
HGA.phi_source = beam_center(2);
HGA.SLx = [0 0 0];
HGA.subarray_geometry = subarray_geometry;
HGA.subarray_window = subarray_window ;

HGA_transmit.theta_source = beam_center(1);
HGA_transmit.phi_source = beam_center(2);
HGA_transmit.SLx = transmit_geometry;


snr = (T * 10^((PSG_transmit+psg_subarray)/10) * Ptransmit *  wavelength^2)./...
    ( (4 * pi)^3 *  10^(NF/10) * kT  * 10^(losses/10) .* (range).^4);
pulse_sv=single((0:(npulses-1))'/(prf));
% count = 0;
% count2 = 1;
% for ii = 1:npulses
%     pulse_sv(ii) = count/prf;
%     if(count2>=transmitters)
%         count = count + 1;
%         count2 = 0;
%     end
%     count2 = count2 + 1;
% end
extra_doppler = zeros(size(pulse_sv));
fc = transmitters/((npulses/transmitters) * 1/prf);
for ii = 1:transmitters
    extra_doppler(ii:transmitters:end) = (transmitters-1-2*(ii)) * fc/2;
end


target  = zeros(length(ranges),npulses,size(geometry ,1));
pulse_replica = pulse_replica/max(max(pulse_replica));
for ii = 1:length(range)
    mask = zeros(size(ranges));
    [junk index] = min(abs(ranges-range));
    index = index(1) - floor(length(pulse_replica));
    if(index<1) 
        index = 1;
    end
    mask(index) = sqrt(snr(ii));

    tf = fftfilt(pulse_replica,mask);

    altitude = (platform_position(3,ii)-target_position(3,ii));
    [EL graz]= sargmtiangles(altitude,sqrt(range(ii).^2 + (altitude * 1000)^2)/1000);
    
    pos = (target_position(1:2,ii)-platform_position(1:2,ii));
    pos = pos/norm(pos);
     uv = [pos(1) * cosd(EL);pos(2) * cosd(EL);sind(EL)];
    doppler_sv =  exp( 2 * pi * 1i * pulse_sv .* ((2 * velocity' * uv)/wavelength  + extra_doppler));

    for ii = 1:transmitters
       
        SV_transmit = 1;%exp(-2 * pi * 1i * mimo_transmit_geometry(:,:,ii)/ wavelength  * [sin(beam_center(1)) * cos(beam_center(2));cos(beam_center(1)) * cos(beam_center(2));sin(beam_center(2))]);
        doppler_sv(:,ii:transmitters:end)  = doppler_sv(:,ii:transmitters:end) *  SV_transmit;
    end
    
    [receive_beam_pattern] = calcBeamPattern(frequency,atan2(pos(1),pos(2)) * 180/pi,EL,[],HGA);
    [transmit_beam_pattern] = calcBeamPattern(frequency,atan2(pos(1),pos(2)) * 180/pi,EL,[],HGA_transmit);
    beam_pattern = receive_beam_pattern .* transmit_beam_pattern;
   
    
    
    
    spatial_sv =    exp(-2 * pi * 1i * geometry * uv / wavelength) * sqrt(beam_pattern);
    C =kron(doppler_sv,spatial_sv.');
    for kk = 1:size(target,1)
        target(kk,:,:) = tf(kk)  * C; 
    end
end




% target_reduced = zeros(length(unambiguous_range),size(target,2),size(target,3));
% for ii = 1:ceil(size(target,1)/length(unambiguous_range))
%     index = (1:length(unambiguous_range)) + (ii-1) * length(unambiguous_range);
%     index(find(index>size(target,1)))= [];
%     index2 = 1:(min([length(ranges) length(index)]));
%     target_reduced(index2,:,:) = target_reduced(index2,:,:)  + target(index,:,:);
% end
% target = target_reduced;