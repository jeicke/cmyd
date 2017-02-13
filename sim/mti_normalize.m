function [D] = mti_normalize(D,method,pulses,shade_doppler)
if(~exist('method','var')||isempty(method))
    method = 'z';
end
channels = size(D,3);
switch lower(method)
    case 'median'
        if(channels~=1)
            for channel = 1:channels
                for l = 1:size(D,2)
                    D(:,l,channel) = D(:,l,channel)./median(abs(D(:,l,channel)));
                end
            end
        else
            for l = 1:size(D,2)
                D(:,l) = D(:,l)./median(abs(D(:,l)));
            end
        end
    case 'fft'
        hf1 = fft([1 -2 1],size(D,2));
        hf1(1) = eps;
        hf = 1./((hf1 ));

        hf = fftshift(hf);
        hf = repmat(hf,size(D,1),[]);
        if(channels~=1)
            for channel = 1:channels
                D(:,:,channel) = D(:,:,channel).*hf;
            end
        else
            D = D.*hf;
        end
    case 'z'
        T=toeplitz([1 -2 1 zeros(1,pulses-3)]',[1 zeros(1,size(D,2)-1)]);
   
        W=dftmtx(size(D,2));
        S = diag([shade_doppler' zeros(1,size(D,2)-length(shade_doppler))]) ;
        for k = 1:size(D,2)
            y(k)=W(:,k)'*S*T'*T*S'*W(:,k);
        end
        
        
        hf = 1./ifftshift(sqrt(abs(y )));
        hf = repmat(hf,size(D,1),[]);
        if(channels~=1)
            for channel = 1:channels
                D(:,:,channel) = D(:,:,channel).*hf;
            end
        else
            D = D.*hf;
        end
end