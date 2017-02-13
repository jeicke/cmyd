%% parameters
Npulses = Npulses * transmitters;
doppler_bins = [-prf/2:2:prf/2];
prf = prf * transmitters;%c/(2 * (clutter_range (2)+5) * 1000);%ceil((350) * 4 * operating_frequency/c);%6 * norm(radar_velocity)/2.4;

v_max = (prf/transmitters) * (c / operating_frequency)/4;
fprintf(1,'***Maximum Unambiguous Velocity %6.2f m/s.  Frequency Resolution %6.3f Hz.\n',v_max,(prf) /Npulses);
%% create geometry of radar

% receive subarray geometry
if(compute_subarray)
    [subarray_geometry]= load_geometry('subarray',[],subarray_elements_width,subarray_elements_height,299792458/5E9);
    % width and height of subarray in meters
    subarray_width =  max(abs(subarray_geometry(:,1))) * 2;
    subarray_height =  max(abs(subarray_geometry(:,3))) * 2;

    % gain of subarray
    PSG_subarray = 10*log10(4 * pi * subarray_width *subarray_height /(299792458/operating_frequency )^2);
    w1 = chebwin(subarray_elements_width,subarray_sidelobes );
    w2 = chebwin(subarray_elements_height,subarray_sidelobes );
else
    subarray_geometry = [0 0 0];
    PSG_subarray = 1;
    subarray_width = 0 ;
    subarray_height = 0;
    w1 = 1;
    w2 = 1;
end
% gain of full array
PSG = 10*log10(4 * pi * width * height /(299792458/operating_frequency )^2);
%subarray weighting

subarray_window = kron(w1/sum(w1),w2/sum(w2));

if(~line_array)
    % now create full receive array geometry
    [fullarray_geometry, w, h]= load_geometry('rectangle',[],width,height,2 * (subarray_width * (1-overlap)));

    %transmit array geometry
    [transmit_geometry, w, h]= load_geometry('rectangle',[],width,height,299792458/operating_frequency);
else
    [fullarray_geometry]= load_geometry('line',[],width,height,299792458/operating_frequency );
    transmit_geometry = fullarray_geometry;
end
if(transmitters>1)

    clear mimo_transmit_geometry;
    oz = 2 * mean(fullarray_geometry(1: size(fullarray_geometry,1)/2,1));
    offset = [(-oz*(transmitters-1)):(oz* 2):(oz* (transmitters-1))];
    mimo_receive_geometry = repmat(fullarray_geometry,transmitters,1);
    for ii = 1:transmitters
        mimo_transmit_geometry(:,:,ii) = [offset(ii)  0 0 ];
        mimo_receive_geometry( (size(fullarray_geometry,1) * (ii-1) + 1): ii * size(fullarray_geometry,1),1) = ...
        mimo_receive_geometry( (size(fullarray_geometry,1) * (ii-1) + 1): ii * size(fullarray_geometry,1),1) + ...
        offset(ii);
    end
    
else
    mimo_transmit_geometry = [];
end

%% derived parameters
unambiguous_range = c/(2 * prf)/1000;

pulse_replica = 1;%m_chirp(chirp_bandwidth,pulse_length,-chirp_bandwidth/2,sampling_rate );
% create targets
resolution = 1;
% range of Azimuths to simulate for clutter
if(~line_array)
    if(compute_radarequation)
        AZ = single( ([-90:resolution:90]) * pi/180);
    else
        AZ = single( ([-10:resolution:10]) * pi/180);
    end
else
     AZ = single( ([-90:resolution:90]) * pi/180);
end
% calculate range spacing in meters
rngspa= 299792458/(sampling_rate)/2;

%ranges from start of recording to horizon
ranges = [clutter_range(1) * 1000:rngspa:(clutter_range(2))*1000];

%set of ranges that occur between pulses
unambiguous_range = [1000 * clutter_range(1):rngspa:(unambiguous_range) * 1000];

%% simulation
fprintf(1,'***Creating Targets\n');
%create target
if(make_targets )
    if(transmitters>1)
        target = target_generate_mimo2(target_position,radar_position,radar_velocity-target_velocity,...
            pulse_replica,ranges , sampling_rate,transmit_beam ,...
            k*T,losses,dwell_time,NF,1,Ptransmit,PSG ,operating_frequency,Npulses,prf,...
            transmit_geometry,fullarray_geometry,subarray_geometry,subarray_window,...
            mimo_transmit_geometry,PSG_subarray ,unambiguous_range,transmitters,range_ambiguous);
    else
        target = target_generate(target_position,radar_position,radar_velocity-target_velocity,...
            pulse_replica,ranges , sampling_rate,transmit_beam ,...
            k*T,losses,dwell_time,NF,1,Ptransmit,PSG ,operating_frequency,Npulses,prf,...
            transmit_geometry,fullarray_geometry,subarray_geometry,subarray_window,PSG_subarray ,unambiguous_range,range_ambiguous);
    end
end
fprintf(1,'***Creating Clutter\n');
% clutter power as a function of angle
   clutter_power = create_power(AZ,pulse_replica,ranges ,sampling_rate,transmit_beam,...
        radar_position(3),k*T,losses,dwell_time,NF,c/operating_frequency,gamma,Ptransmit,...
        PSG,operating_frequency,fullarray_geometry,transmit_geometry,subarray_geometry,subarray_window,PSG_subarray ,compute_radarequation,power_clipping);
%doppler for each clutter ring
if(ISBISTATIC)
    [dopplers, EL]= create_doppler_bi(AZ,radar_velocity, ranges,operating_frequency,radar_position(3));
    if(~ALLOWWALK)
        %get doppler at peak illumination
        HGA_transmit.theta_source =transmit_beam(1);
        HGA_transmit.phi_source = transmit_beam(2);
        HGA_transmit.SLx = transmit_geometry;
        [theory_matrix_transmit] = calcBeamPattern(operating_frequency,AZ * 180/pi,transmit_beam(2),[],HGA_transmit);
        AZMAX = 0;
        [doppleraz0, EL] = create_doppler_bi(AZMAX,radar_velocity, ranges,operating_frequency,radar_position(3));
    else
        doppleraz0 = [];
    end
else
    [dopplers, EL]= create_doppler(AZ,radar_velocity, ranges,operating_frequency,radar_position(3));
    if(~ALLOWWALK)
    end
end


% integrate clutter patches around circle
if(make_clutter)
    if(isempty(range_gate))
        if(transmitters>1)
            [clutter]= create_clutter_mimo2(AZ,clutter_power,dopplers,(operating_frequency)/299792458,fullarray_geometry,mimo_transmit_geometry,EL,Npulses,prf,ranges,unambiguous_range, power_threshold, transmitters,transmit_beam ,range_ambiguous );
        else
            [clutter]= create_clutter4(AZ,clutter_power,dopplers,(operating_frequency)/299792458,fullarray_geometry,EL,Npulses,prf,ranges,unambiguous_range, power_threshold,range_ambiguous);
        end
    else
         
        [R2, Rn]= create_clutter_r(AZ,clutter_power,dopplers,(operating_frequency)/299792458,fullarray_geometry,mimo_transmit_geometry,EL,Npulses,prf,ranges,unambiguous_range, power_threshold, transmitters,transmit_beam ,range_ambiguous,range_gate,dopplers);
        %          clear R;
        %          R = zeros(size(s,1),size(s,1),size(s,3));
        %          for ii = 1:length(S)
        %              s = S{ii};
        %              R(:,:,ii) = s * s';
        %          end
        if(transmitters<2)



            [V] = steering_vector(fullarray_geometry, receive_beam(1) * pi/180, ...
                receive_beam(2) * pi/180, (operating_frequency)/299792458 );
            clear R;
            for jj =1:size(R2,3)

                R(:,:,jj) = R2(:,:,jj) * R2(:,:,jj)';
            end
            R = R/size(R2,2);
            if(compute_stap)
                R3 = reshape(R2(:),Npulses,size(fullarray_geometry,1),[]);
                R3 = permute(R3,[3,1,2]);
                clear R;
                for jj =1:size(Rn,3)

                    R(:,:,jj) = Rn(:,:,jj) * Rn(:,:,jj)';
                end
                R = R/size(R2,2);
                %                  clear R;
                %                  R2 = reshape(R4(:),Npulses*size(fullarray_geometry,1),[]);
                %                 for jj =1:size(R2,3)
                %
                %                     R(:,:,jj) = R2(:,:,jj) * R2(:,:,jj)';
                %                 end
                %                 R = R/size(R2,2);
                sinr = sinr_predoppler(R3,Npulses,size(fullarray_geometry,1),doppler_bins,V,pre_doppler_response',prf,R);
            elseif(compute_post_doppler_stap)
                R3 = reshape(R2(:),Npulses,size(fullarray_geometry,1),[]);
              
                sinr = sinrR_postdoppler(R3,Npulses,size(fullarray_geometry,1),doppler_bins,V,prf);
            else

                sinr = sinrR(R,Npulses,size(fullarray_geometry,1),doppler_bins,V,prf);
            end
        else

            [V] = steering_vector(mimo_receive_geometry, receive_beam(1) * pi/180, ...
                receive_beam(2) * pi/180, (operating_frequency)/299792458 );
            if(compute_stap)
                clear R;
                clear D2;
                for kk =1:size(R2,3)
                    R3 = reshape(R2(:,:,kk),Npulses,size(fullarray_geometry,1),[]);
                    R3 = permute(R3,[3,1,2]);
                    [D] = doppler_filter(R3,Npulses,'none',70,false)/transmitters;

                    clear D4;
                    for ii = 1:transmitters
                        jj = transmitters-(ii-1);
                        D4(:,:,(ii-1) * size(fullarray_geometry,1)+1:ii*size(fullarray_geometry,1)) = D(:,((Npulses/transmitters) * (jj-1)  + 1):((Npulses/transmitters) * jj),:);
                    end
                   
                    [D] = conj(inverse_doppler_filter(D4));
                    
                    D2(:,:,:,kk) = D;
                    S = reshape(D,size(D4,1),Npulses/transmitters,size(mimo_receive_geometry,1));
                    S = permute(S,[1, 3, 2]);
                    S = reshape(S(:),size(D4,1),[]);

                    R(:,:,kk) = S'*S;

                end
                sinr = sinr_predoppler(conj(D2),Npulses/transmitters,size(mimo_receive_geometry,1),doppler_bins,V,pre_doppler_response',prf/transmitters,R);
            elseif(compute_post_doppler_stap)

            else
                clear R;
                for kk =1:size(R2,3)
                    R3 = reshape(R2(:,:,kk),Npulses,size(fullarray_geometry,1),[]);
                    R3 = permute(R3,[3,1,2]);
                    [D] = doppler_filter(R3,Npulses,'none',70,false)/transmitters;
                    clear D4;
                    for ii = 1:transmitters
                        jj = transmitters-(ii-1);
                        D4(:,:,(ii-1) * size(fullarray_geometry,1)+1:ii*size(fullarray_geometry,1)) = D(:,((Npulses/transmitters) * (jj-1)  + 1):((Npulses/transmitters) * jj),:);
                    end
                
                    [D] = conj(inverse_doppler_filter(D4));
                    S = reshape(D,size(D4,1),Npulses/transmitters*size(mimo_receive_geometry,1),[]);
                    R(:,:,kk) = S'*S;

                end
                sinr = sinrR(R,Npulses/transmitters,size(mimo_receive_geometry,1),doppler_bins,V,prf/transmitters);
            end
            a = 1;
        end
         return;
    end
end

% add noise
fprintf(1,'***Adding noise \n');
noise = (1/sqrt(2))*(randn(size(clutter)) + 1i * randn(size(clutter)));
% compute receive beam
D = zeros(size(clutter));
if(make_targets) 
    D = D + target; 
end
if(make_clutter) 
    D = D + clutter;
end
if(make_noise) 
    if(transmitters>1)
        f = Npulses/(transmitters * transmitters);
    else
        f = Npulses;
    end
    D = D + sqrt(pulse_length * chirp_bandwidth) *  sqrt(f) * noise;
end

%%
fprintf(1,'***Processing radar\n');
[Ds ] = pulse_compression(D,pulse_replica,'cheb',70);
Ds = Ds(length(pulse_replica):end,:,:);
[D] = doppler_filter(Ds,Npulses,'cheb',70,false)/transmitters;
if(transmitters>1)
    clear D4;
    
    for ii = 1:transmitters
        jj = transmitters-(ii-1);
         D4(:,:,(ii-1) * size(fullarray_geometry,1)+1:ii*size(fullarray_geometry,1)) = D(:,((Npulses/transmitters) * (jj-1)  + 1):((Npulses/transmitters) * jj),:);
    end


    [V] = steering_vector(mimo_receive_geometry, receive_beam(1) * pi/180, ...
       receive_beam(2) * pi/180, (operating_frequency)/299792458 );
    [V_sweep] = steering_vector(mimo_receive_geometry, az * pi/180, ...
       receive_beam(2) * ones(size(az)) * pi/180, (operating_frequency)/299792458 );
   [V_sweep_el] = steering_vector(mimo_receive_geometry, receive_beam(1) * pi/180, ...
       el * pi/180, (operating_frequency)/299792458 );
    D = D4;
    D_mimo = D;
else
   
    [V] = steering_vector(fullarray_geometry, receive_beam(1) * pi/180, ...
       receive_beam(2) * pi/180, (operating_frequency)/299792458 );
    [V_sweep] = steering_vector(fullarray_geometry, az * pi/180, ...
       receive_beam(2) * ones(size(az)) * pi/180, (operating_frequency)/299792458 );
     [V_sweep_el] = steering_vector(fullarray_geometry, receive_beam(1) * pi/180, ...
       el * pi/180, (operating_frequency)/299792458 );
   D_simo = D;
end
[Dbf ] = beamform(D,1:size(D,3),V'.');
[Dbf_sweep ] = beamform(D,1:size(D,3),V_sweep'.');
[Dbf_sweep_el ] = beamform(D,1:size(D,3),V_sweep_el'.');

f = prf/transmitters * (-size(Dbf,2)/2:size(Dbf,2)/2-1)/(size(Dbf,2)/2);
imagesc(fftshift((20*log10(abs(Dbf(:,:,1)))),2))
caxis([0 90])
colorbar
drawnow;

if(transmitters>1)
    Dbf_mimo = Dbf;
    Dbf_sweep_mimo = Dbf_sweep;
    Dbf_sweep_el_mimo = Dbf_sweep_el;
    f_mimo = f;
else
    Dbf_simo = Dbf;
    Dbf_sweep_simo = Dbf_sweep;
    Dbf_sweep_el_simo = Dbf_sweep_el;
    
    
    f_simo = f;
end

if(compute_stap)
    
    fprintf(1,'***Processing radar STAP\n');
    if(transmitters>1)
        [D] = doppler_filter(Ds,Npulses,'none',70,false);
        for ii = 1:transmitters
            jj = transmitters-(ii-1);
            D4(:,:,(ii-1) * size(fullarray_geometry,1)+1:ii*size(fullarray_geometry,1)) = D(:,((Npulses/transmitters) * (jj-1)  + 1):((Npulses/transmitters) * jj),:);
        end
        [Ds] = inverse_doppler_filter(D4);
        sensors = size(mimo_receive_geometry,1);
    else
        sensors = size(fullarray_geometry,1);
    end
    [Dbf2 w] = beamform_stap(Ds,V'.',[0 1 0]',20,10,1E-4,20);
    [x ii] = min(abs(receive_beam(2)- EL * 180/pi));
    
    [sinr sinr2] = calculate_sinr(w,  Ds,V'.',size(Ds,2),sensors ,prf/transmitters ,20,10,20,ii:ii+1);
    

    [Dbf2 shade_doppler] = doppler_filter(Dbf2,Npulses,'cheb',70);
    [Dbf2] = mti_normalize(Dbf2,'median');
    figure;
    p = 20 * log10(abs(Dbf2));
    imagesc(f, unambiguous_range/1000,p)
    caxis([0 70])
    colorbar
    drawnow;
    if(transmitters>1)
        Dbf_stap_mimo = Dbf2;
        sinr_stap_mimo = sinr2;
    else
        Dbf_stap_simo = Dbf2;
        sinr_stap_simo = sinr2;
    end
end
if(compute_post_doppler_stap)
    
    fprintf(1,'***Processing radar STAP\n');
    if(transmitters>1)
        [D] = doppler_filter(Ds,Npulses,'none',70,false);
        for ii = 1:transmitters
            jj = transmitters-(ii-1);
            D4(:,:,(ii-1) * size(fullarray_geometry,1)+1:ii*size(fullarray_geometry,1)) = D(:,((Npulses/transmitters) * (jj-1)  + 1):((Npulses/transmitters) * jj),:);
        end
        Ds = D4;
        sensors = size(mimo_receive_geometry,1);
    else
        sensors = size(fullarray_geometry,1);
        Ds = doppler_filter(Ds,Npulses,'cheb',70,false);
    end
    [dbf2 w sinr] = beamform_postdoppler_stap(Ds,V'.',20,10,1E-4,20);
   
    [x ii] = min(abs(receive_beam(2)- EL * 180/pi));
    
    
    

    
    figure;
    p = 20 * log10(abs(dbf2));
    imagesc(f, unambiguous_range/1000,p)
    caxis([0 70])
    colorbar
    drawnow;
    if(transmitters>1)
        Dbf_stap_mimo = dbf2;
        sinr_stap_mimo = sinr;
    else
        Dbf_stap_simo = dbf2;
        sinr_stap_simo = sinr;
    end
end

clear clutter;
clear target;
clear noise;