function displayparameters( Parameters )
%PRINTPARAMETRS Prints the parameters used in simulation to screen
fprintf(1,'\tBaseband sampling rate %d MHz\n',Parameters.basebandSamplingRate/1E6);
fprintf(1,'\tCarrier Frequency rate %d MHz\n',Parameters.carrierFrequency/1E6);
fprintf(1,'\tOut files will start with the prefix %s\n',Parameters.fileNameBase  );
fprintf(1,'\tTotal time in simulation is %3.3f seconds\n',Parameters.time(end));
fprintf(1,'\tMax time delay in simulation is %3.5f microseconds\n',Parameters.maxSystemDelay*1E6);
fprintf(1,'\tMinimum time delay in simulation is %3.5f microseconds\n',Parameters.minSystemDelay*1E6);
fprintf(1,'\tTotal objects in simulation is %d\n',size(Parameters.track,3));
if(Parameters.rangeCorrect)
    fprintf(1,'\tComputing range attenuation\n');
else
    fprintf(1,'\tNot computing range attenuation\n');
end
end

