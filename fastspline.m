function [timeseries state ] = fastspline(state,Fs, t,y,cast)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
timeseries = [];
if(state.init)
    h = 1/Fs;
    % create interpolation coeffs
    state.t2 = [t t(end) + h];
    y2 = [y y(end)];
    n = length(state.t2)-2;
    D = sparse(1:n,1:n,4*ones(1,n),n,n);
    E = sparse(2:n,1:n-1,ones(1,n-1),n,n);
    S = E+D+E';
    m = length(state.t2);
    rhs = sparse(6/h^2 * (y2(1:m-2)-2*y2(2:m-1) + y2(3:m)));
    x = full(S\rhs.');
    M = [0; x; 0].';
    state.a = cast((M(2:end)-M(1:end-1))/(6*h));
    state.b = cast(M(1:end-1)/2);
    state.c = cast((y2(2:end)-y2(1:end-1))/h-(M(2:end)+2 * M(1:end-1)) * h/6);
    state.d = cast(y2(1:end-1));
    state.t2 = cast(state.t2);
    state.init = false;
else
     h = 1/Fs;
    [v i_start] = min(double(abs(t(1) - state.t2 )));
    
    index = floor(t * Fs) -floor(t(1) * Fs);
    index = ( i_start+index);
    s =  t-state.t2(index(1));
    if(s(1)<0)
        index  = index -1;
        s =  t-state.t2(index);
    elseif(s(1)>1)
        index  = index +1;
        s =  t-state.t2(index);
    else
       
        s =  t-state.t2(index);
    end

    index = index -1;
    timeseries = ((state.a(index).*s+state.b(index)).*s+state.c(index)).*s+state.d(index);

end

end

