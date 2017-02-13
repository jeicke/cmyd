clear Pxx
for angle = [1:3]
    %%
    filebase = sprintf('//bigsazenas/data/sprdn/statictank_%d',angle);
    load(sprintf('//bigsazenas/data/sprdn/statictank_%d_Parameters.mat',angle));
    pulse_replica=m_chirp(Parameters.timeSeriesParameters(1),Parameters.timeSeriesParameters(2),0,Parameters.basebandSamplingRate);
    
    [systemDelay(2), systemDelay(1)] = finddelays( Parameters.track,[0 0 0],Parameters.C,Parameters.timeSeriesParameters(2) );
    %[P, blockTimes] =processradarblock( filebase,pulse_replica,systemDelay);
    Fs = Parameters.basebandSamplingRate;
    C = Parameters.C;
    Z = floor((systemDelay(1))*Fs);
    Tz = Z/Fs;
    ranges = ([1:1:Parameters.blockSize]-1) * C/Fs * .5+Tz * C * .5-C/Fs*length(pulse_replica)/4+C/Fs/2;
    
    
    
    
    %% Point spread function
    trackCenter = [ 0 0 0];
    xi = [-15:.05:15] + trackCenter(1);
    yi = [-15:.05:15] +  trackCenter(2);
    zi = [0 ] +  trackCenter(3);
    track = Parameters.track;
    prf = Parameters.timeSeriesParameters(3);
    [Px] =bpm( filebase,prf,pulse_replica,ranges,track,[0 0 0],xi,yi,zi);
    %%
    save(sprintf('//bigsazenas/data/sprdn/statictank_%d_sarimage.mat',angle),'Px','xi','yi','zi','prf','pulse_replica','track','Parameters','-v7.3' );
    figure
    colormap jet
    imagesc(-xi/1000,yi,20*log10(abs(squeeze(Px))));
    nf = median(20*log10(abs(Px(:))));
    caxis([nf-5 max(20*log10(abs(squeeze(Px(:)))))]);
    ylabel('CROSS-RANGE (M)')
    xlabel('RANGE (KM)');
    hold on;
    colorbar
    plot(squeeze(track(1,1,3:end)),squeeze(track(2,1,3:end)),'ro')
    Pxx{angle} = Px;
end