function [s state] = chirppulse( times,state, parameters,Fs )
%tone produces a sine wave at a given frequency
%  parameters is frequency if unspecified defaults to 25 MHz

if(~exist('state')|| ~isfield(state,'savedTs'))

    %sin( 2 * pi * 1 * times);%
      state.frequency = parameters(3);
    state.chirp=m_chirp(parameters(1),parameters(2),0,Fs);
    if(length(state.chirp>6))
     %   state.chirp=state.chirp.*hann(length(state.chirp));
    end
    pulse_length = round(Fs/ state.frequency);
    state.chirp = [state.chirp.' zeros(1,pulse_length-length(state.chirp))];
    state.savedTs = zeros(size(times));%sin(2 * pi * state.frequency * times);%sin( 2 * pi * 1 * times);%
   
    state.analytic = false;
    state.times = times-times(end)-1/Fs;
 
end

if(state.create)
    index = find( (state.times(end)- abs(state.times))> state.maxSystemDelay);
    state.times(index) = [];
    state.savedTs(index) = [];

    s1 = round(times * Fs);%sin(2 * pi * state.frequency * times);%
    s1 = state.chirp(mod(s1,length(state.chirp))+1);% double((s1/state.frequency)==floor((s1/state.frequency)));
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


