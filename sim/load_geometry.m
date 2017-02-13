function [geometry w h]= load_geometry(array,channels,width,height,wavelength)
geometry = [];
switch lower(array)
    case {'line'}
         channels= ceil(width/(wavelength))/2;
        for ii = 1:channels
            geometry(ii,:) = [width-(ii-1)* wavelength 0 0];
        end
        geometry(:,1) = geometry(:,1) - mean(geometry(:,1));
    case {'subarray'}
        w = width;
        h = height;
        count  = 1;
        for ii = 1:w
            for jj = 1:h
                 geometry(count,:) =  [width-(ii-1)* wavelength/2 0 height-(jj-1)* wavelength/2];
                 count  = count + 1;
            end
        end
        geometry(:,1) = geometry(:,1) - mean(geometry(:,1));
        geometry(:,3) = geometry(:,3) - mean(geometry(:,3));
    case {'rectangle'}
        w = ceil(width/(wavelength/2))/2;
        h = ceil(height/(wavelength/2))/2;
        count  = 1;
        for ii = 1:w
            for jj = 1:h
                 geometry(count,:) =  [width-(ii-1)* wavelength 0 height-(jj-1)* wavelength];
                 count  = count + 1;
            end
        end
        geometry(:,1) = geometry(:,1) - mean(geometry(:,1));
        geometry(:,3) = geometry(:,3) - mean(geometry(:,3));
    case 'uesa'
         [Element_Angle Element_Radii ] = Radar_Globals();
         geometry = [sin(Element_Angle) .* Element_Radii; cos(Element_Angle) .* Element_Radii; zeros(size(Element_Angle))];
         % now add aux array geometry
         channel_mapping = [1 2 4 5:14 16:26];
         geometry(3,28:54) = geometry(3,28:54) + 10*12 * .0254;
         for ii = 1:26
             geometry(:,ii) = [0;0;0];
         end
         
         
         geometry(:,channel_mapping(1:12)) = 12 * [zeros(1,12);1:12;zeros(1,12)] * .0254;
         geometry(2,channel_mapping(1:12)) = geometry(2,channel_mapping(1:12)) - mean( geometry(2,channel_mapping(1:12)),2);
         geometry(3,channel_mapping(1:12)) = 12* .0254;
         geometry(:,channel_mapping(13:24)) = 12 * [zeros(1,12);1:12;zeros(1,12)] * .0254;
         geometry(2,channel_mapping(13:24)) = geometry(2,channel_mapping(13:24)) - mean( geometry(2,channel_mapping(13:24)),2);
         geometry(3,channel_mapping(13:24)) = -12* .0254;
         geometry(3,channel_mapping) = geometry(3,channel_mapping) - 10*12 * .0254;
         geometry = geometry';
         
%         [Element_Angle Element_Radii ] = Radar_Globals();
%          geometry = [sin(Element_Angle) .* Element_Radii; cos(Element_Angle) .* Element_Radii; zeros(size(Element_Angle))];
%          % now add aux array geometry
%          channel_mapping = [1 2 4 5:14 16:26];
%          for ii = 1:26
%              geometry(:,ii) = [0;0;0];
%          end
%          
%          
%          geometry(:,channel_mapping(1:12)) = 12 * [zeros(1,12);1:12;zeros(1,12)] * .0254;
%          geometry(2,channel_mapping(1:12)) = geometry(2,channel_mapping(1:12)) - mean( geometry(2,channel_mapping(1:12)),2);
%          geometry(3,channel_mapping(1:12)) = 12* .0254;
%          geometry(:,channel_mapping(13:24)) = 12 * [zeros(1,12);1:12;zeros(1,12)] * .0254;
%          geometry(2,channel_mapping(13:24)) = geometry(2,channel_mapping(13:24)) - mean( geometry(2,channel_mapping(13:24)),2);
%          geometry(3,channel_mapping(13:24)) = -12* .0254;
%          geometry = geometry';
    case 'combined'
        [Element_Angle Element_Radii ] = Radar_Globals();
         geometry = [sin(Element_Angle) .* Element_Radii; cos(Element_Angle) .* Element_Radii; zeros(size(Element_Angle))];
         % now add aux array geometry
         channel_mapping = [1 2 4 5:14 16:26];
         geometry(3,28:54) = geometry(3,28:54) + 10*12 * .0254;
         for ii = 1:26
             geometry(:,ii) = [0;0;0];
         end
         
         
         geometry(:,channel_mapping(1:12)) = 12 * [zeros(1,12);1:12;zeros(1,12)] * .0254;
         geometry(2,channel_mapping(1:12)) = geometry(2,channel_mapping(1:12)) - mean( geometry(2,channel_mapping(1:12)),2);
         geometry(3,channel_mapping(1:12)) = 12* .0254;
         geometry(:,channel_mapping(13:24)) = 12 * [zeros(1,12);1:12;zeros(1,12)] * .0254;
         geometry(2,channel_mapping(13:24)) = geometry(2,channel_mapping(13:24)) - mean( geometry(2,channel_mapping(13:24)),2);
         geometry(3,channel_mapping(13:24)) = -12* .0254;
         geometry(3,channel_mapping) = geometry(3,channel_mapping) - 10*12 * .0254;
         geometry = geometry';
    otherwise
        disp('Unknown method.')
end