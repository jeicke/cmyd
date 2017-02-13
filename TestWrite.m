function  [fp] = TestWrite(timeSeries,firstFlag,filename,Fs,Fc,channels,fp)
MAX_SHORT = 32768;
MAX_V = 2*sqrt(2);

%timeseries is samples X channels
%check and see if this is first call and open file
try
    if(firstFlag)
        fp = fopen(filename, 'w');
        fwrite(fp,Fs,'float32');
        fwrite(fp,Fc,'float32');
        fwrite(fp,channels,'float32');
    else
        if(~exist('fp','var')||isempty(fp))
            fp = fopen(filename, 'a');
        end
    end
catch
    rethrow(lasterror)
end
z = [real(timeSeries); imag(timeSeries)];
z = z(:);
fwrite(fp,z,'float32');


%fclose(fp);


%end of function%
