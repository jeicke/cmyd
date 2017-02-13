function [data Fs Fc channels] = testread(filename,blocks_offset,blocks_read,blockLength)
fp = fopen(sprintf('%s_%s',filename,'received.bin'), 'r');
% read header
temp = fread(fp,3,'float32');
Fs = temp(1);
Fc = temp(2);
channels = temp(3);

if(~exist('blocks_offset','var')||isempty(blocks_offset))
    blocks_offset = 0;
    blocks_read = inf;
end
if(blocks_offset == -1)
    fseek(fp,0,'eof');
    pos = ftell(fp);
    data = (pos-4 * 3)/8/channels;
else
    fseek(fp,blocks_offset * channels + 12,'bof');
    temp = fread(fp,blocks_read*2* channels,'float32=>single');
    data = temp(1:2:end) + 1i * temp(2:2:end);
    
%     temp = fread(fp,[2 blocks_read* channels],'float32=>single');
%     data = temp(:,1) + 1i * temp(:,2);
    
    data = reshape(data,blockLength,channels,[]);
end
fclose(fp);