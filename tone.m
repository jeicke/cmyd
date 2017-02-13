function [s state] =tone( times,state, parameters,Fs )
%tone produces a sine wave at a given frequency
%  parameters is frequency if unspecified defaults to 25 MHz

if(~exist('state')|| ~isfield(state,'savedTs'))

    %sin( 2 * pi * 1 * times);%
    state.frequency = parameters;
    state.savedTs = zeros(size(times));%sin(2 * pi * state.frequency * times);%sin( 2 * pi * 1 * times);%
   
    state.analytic = false;
    state.times = times-times(end)-1/Fs;
 
end

if(state.create)
    index = find( (state.times(end)- abs(state.times))> state.maxSystemDelay);
    state.times(index) = [];
    state.savedTs(index) = [];

    s1 = cos(2 * pi * state.frequency * times);%
    if(state.frequency > 100000)
        s_new = sqrt(length(s1))/(norm(s1)) *  s1;
        
    else
         s_new = s1;
    end
    state.savedTs = [state.savedTs s_new];
    state.times = [state.times times];
    
    s = s_new;
else

    [v i_start] = min(abs(times(1) - state.times ));
    s = state.savedTs(i_start:1:length(times) + i_start-1);

end
end

%pulse_replica = 1;%m_chirp(chirp_bandwidth,pulse_length,-chirp_bandwidth/2,sampling_rate ).';
%pulse_replica = pulse_replica .*chebwin((length(pulse_replica)),70).';

% if(~exist('state')|| ~isfield(state,'frequency'))
%     if(exist('parameters','var')&&~isempty(parameters))
%         state.frequency = parameters;
%         state.analytic = true;
%     end
% end
% s = sin(2 * pi * state.frequency * times).';



%tone produces a sine wave at a given frequency
%   if frequency is unspecified defaults to 25 MHz