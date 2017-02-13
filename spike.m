function [s state] =spike( times,state, parameters,Fs )
%tone produces a sine wave at a given frequency
%  parameters is frequency if unspecified defaults to 25 MHz

if(~exist('state')|| ~isfield(state,'savedTs'))

    %sin( 2 * pi * 1 * times);%
    state.frequency = parameters(3);
    state.spike = 1;
    pulse_length = round(Fs/ state.frequency);
    state.spike = [1 zeros(1,pulse_length-1)];
    state.savedTs = zeros(size(times));%sin(2 * pi * state.frequency * times);%sin( 2 * pi * 1 * times);%
   
    state.analytic = false;
    state.times = times-times(end)-1/Fs;
 
end

if(state.create)
    index = find( (state.times(end)- abs(state.times))> state.maxSystemDelay);
    state.times(index) = [];
    state.savedTs(index) = [];

    s1 = round(times * Fs);%sin(2 * pi * state.frequency * times);%
    s1 = state.spike(mod(s1,length(state.spike))+1);
%     if(state.frequency > 100000)
%         s_new = sqrt(length(s1))/(norm(s1)) *  s1;
%         
%     else
%          s_new = s1;
%     end
    s_new = s1;
    state.savedTs = [state.savedTs s_new];
    state.times = [state.times times];
    
    s = s_new;
else

    [v i_start] = min(abs(times(1) - state.times ));
    s = state.savedTs(i_start:1:length(times) + i_start-1);

end
end



% if(~exist('state')|| ~isfield(state,'frequency'))
%     if(exist('parameters','var')&&~isempty(parameters))
%         state.frequency = parameters;
%         state.analytic = true;
%     end
% end
% s = sin(2 * pi * state.frequency * times).';



%tone produces a sine wave at a given frequency
%   if frequency is unspecified defaults to 25 MHz