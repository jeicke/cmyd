radsim(filebase,time,track,0,...
       'timeseries' ,@chirppulse,'parameters',[pulseBandwidth pulseLength prf],...
       'fs',Fs,'fc',Fc,'transmitgeometry',geometry,...
       'transmitorientation',orientation,'clutter',reflextivity',...
       'transmitshading',w); 