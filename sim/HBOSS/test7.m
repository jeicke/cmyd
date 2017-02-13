width = [2.6 7.3 .5 .6];
height = [.4 .75 .2 .6];
power = [1610 300 165 60];
altitude = [20 20 5 5 ];
beams = [32 85 6 7];
range = [25 300;25 300;7 100;7 100];
ii = 2;
 f = [5.5 10 17] * 1E9;
        cpi = [3 5 7];
        snr = [13.5 11.5 10.5];
        elev_ranges{1} = [ 300.0000   30.1650   25];
            elev_ranges{2} = [300.0000   47.4227   26.1997   25];
            elev_ranges{3} = [300.0000   74.0697   43.2877   30.6287   25];
for jj = 1:length(f)
      
         [doppler_dwell_time_beam power_dwell_time az_beams power_dwell_time_band doppler_dwell_time_band X Y] = gmtiscript(f(jj),snr(jj),cpi(jj),power(ii),altitude(ii),range(ii,:),width(ii),height(ii),...
             [],[],[],[],[],[],[],[],elev_ranges{jj},[],true);
          ratio2 = ceil(doppler_dwell_time_band./power_dwell_time_band);
        [I] = find(ratio2(:)>84);
        ratio2 = ratio2(:);
        ratio2(I) = 84;
        ratio2 = reshape(ratio2,276,[]);
        figure;
        
        pcolor(X,Y,ratio2)
        shading flat;
        axis equal
        axis tight
        xlabel('Y KM')
        ylabel('X KM')
        h = colorbar;
        ylabel(h,'AZIMUTHAL BEAMS');
        if(jj == 1)
            label = sprintf('C:/figures/fan_jstar_elaz_5_ghz');
        elseif(jj==2)
            label = sprintf('C:/figures/fan_jstar_elaz_10_ghz');
        else
            label = sprintf('C:/figures/fan_jstar_elaz_17_ghz');
        end
        
        print( gcf, '-djpeg', label)
end                                                       