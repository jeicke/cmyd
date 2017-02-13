if(1)
    tic 
    cast = @GPUsingle;
    
   % cast = @double;
    %D2 = single(D2);
    
   % cast = @single;
    
    N = size(D2,1);

    M = 1:size(D2,2);
    pulse_replica=m_chirp(pulseBandwidth,pulseLength,0,Fs);
    weighting  =  chebwin(length(pulse_replica),50) ;
    weighting =double(weighting);
    pulse_replica =  double(pulse_replica./norm(pulse_replica));
    pulse_replica = pulse_replica .*weighting;
    P = single(fft(pulse_replica ,N).');
    P = cast(repmat(P,length(M),[])).';
    
    f = single(ifftshift((-N/2:(N/2-1))'/N*Fs));
    K = cast(-1i * 2 * pi  * (f - Fc));
    
    
    LLL = [];
    
    
    [ trackt] =single(maketrack(time,time/10,transmitterVelocity,receiverVelocity ,transmitterPosition,receiverPosition ));
    for ii = 1:size(trialpoint,2);
        %scatteredpoint = [rand(1)-.5 rand(1)-.5 0] * radius ; %relative to transmitetr
        trackt(:,:,ii+2) = trackt(:,:,1);
       %  trackt(1:3,:,ii+2) = 0;
        trackt(1,:,ii+2) = trackt(1,:,ii+2) + trialpoint(1,ii);
        trackt(2,:,ii+2) = trackt(2,:,ii+2) + trialpoint(2,ii);
        trackt(3,:,ii+2) = trackt(3,:,ii+2) + trialpoint(3,ii);
    end
    
    sar = (single(zeros(1,size(trialpoint,2),size(D2,3))));

   
    blockTimes2 = single(blockTimes2);
    ii = M;
    CC = cast(ones(1,size(D2,1)));
    BB = cast(ones(length(M),1));
   
%    setComplex(BB);

   % setComplex(CC);

    for kk = 1:max([1 (size(D2,3)-1)])
       
            Dm = fft(squeeze(D2(:,ii,kk)),N);
            Dm2 = cast(Dm);%cast(repmat(Dm.',size(trialpoint,2),[]).');
            Dm2 = Dm2(end:-1:1,:);
            %Dm2 = ifft(Dm2);
            Dm2 = P.*Dm2;

            count = 1;
            l = 1;
            M = cast(zeros(size(Dm2)));
           
           for jj = 1:size(trialpoint,2)

                
                %-.75*C * length(pulse_replica)/Fs )+C*1/Fs;
                
                [r_transmitter_scatter]= (computeRange(trackt(:,:,[1 2 jj+2]),  blockTimes2(ii,kk)',1));
                
                % next comput range from scatterer to receiver and the angles
                %bsxfun(@plus,scatterPos,-antPos)
                [r_scatter_ant]=(computeRange(trackt(:,:,[1 2 jj+2]), blockTimes2(ii,kk)',2));
                
                
              
                
                % add path from transmitter to scatter and from scatter to receiver  together for scattered path
                
                rangelook = cast((r_transmitter_scatter(2:end,:) + r_scatter_ant(2:end,:))-1*C/Fs+rangeShift)/C;% * length(pulse_replica)/Fs;
               
        

                sar(:,jj,kk) = sar(:,jj,kk) +  single( CC*((Dm2.*exp(K *rangelook))*BB) );
               count = count + 1;
               if(count>100)
                   fprintf('Percent done %2.2f\n',100*min([1 jj/size(trialpoint,2)]))
                   count = 1;
               end
               % sar(:,jj) = sar(:,jj) +   (DD.*Dm2);
            end
       
        %MMM = reshape(single(sar),51,[]);
    end
    toc
end
