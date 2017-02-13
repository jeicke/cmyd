function clutter_power = create_power(AZ,pulse_replica,range,BW,beam_center,...
                                      altitude,kT,losses,T,NF,wavelength,gamma,Ptransmit,PSG_transmit,frequency,...
                                      fullarray_geometry,transmit_geometry,subarray_geometry,subarray_window,psg_subarray,compute_radarequation,power_clipping)
% calculate ground patch area



clear HGA;
HGA.theta_source =beam_center(1);
HGA.phi_source = beam_center(2);
HGA.SLx = [0 0 0];
HGA.subarray_geometry = subarray_geometry ;
HGA.subarray_window = subarray_window ;
resolution = 1;
[theory_matrix_receive] = calcBeamPattern(frequency,AZ * 180/pi,beam_center(2),[],HGA);

HGA_transmit.theta_source =beam_center(1);
HGA_transmit.phi_source = beam_center(2);
HGA_transmit.SLx = transmit_geometry;
[theory_matrix_transmit] = calcBeamPattern(frequency,AZ * 180/pi,beam_center(2),[],HGA_transmit);

theory_matrix = theory_matrix_receive .* theory_matrix_transmit;
rngspa= 299792458/(BW)/2;


clutter_ranges = range;
[EL graz]= sargmtiangles(altitude,sqrt(clutter_ranges.^2 + (altitude * 1000)^2)/1000);
HGA.theta_source =beam_center(1);
HGA.phi_source = beam_center(2);
HGA.SLx = [0 0 0];
HGA.subarray_geometry = subarray_geometry ;
HGA.subarray_window = subarray_window ;
[theory_matrix_EL_receive az de HGA] = calcBeamPattern(frequency,beam_center(1),EL,[],HGA);
[theory_matrix_EL_transmit] = calcBeamPattern(frequency,beam_center(1),EL,[],HGA_transmit);
theory_matrix_EL = theory_matrix_EL_receive.* theory_matrix_EL_transmit;

% index = find(clutter_ranges <clutter_range(1)*1000|
% clutter_ranges>clutter_range(2)*1000);plot(
% index = index - floor(length(pulse_replica));
% index(find(index<1)) = [];
mask = theory_matrix_EL';
%mask(index) = 0;
graz = abs(graz);
clutterarea =  sqrt(2 * clutter_ranges.^2 * (1 - cosd(resolution))) * rngspa ./ cosd(graz) ;%range*rngspa./cos(graz*pi/180);
clutterarea = clutterarea .* sind(graz).*(10.^(gamma/10));
snr = (T * 10^((PSG_transmit+psg_subarray)/10) * Ptransmit *  wavelength^2)./...
    ( (4 * pi)^3 *  10^(NF/10) * kT  * 10^(losses/10) .* (clutter_ranges).^4);


snr = snr .* clutterarea;

snr(find(10*log10(snr)>power_clipping)) = 10^(power_clipping/10);
mask2 = [zeros(1,length(pulse_replica)) mask zeros(1,length(pulse_replica))];
tf = fftfilt(double(pulse_replica),double(mask2));
tf = tf(length(pulse_replica)+1:(end-length(pulse_replica)));
if(compute_radarequation)
    clutter_power = kron(tf.*snr,theory_matrix.');
else
    clutter_power = kron(tf.*max(snr),ones(size(theory_matrix.')));
end
a = 1;
