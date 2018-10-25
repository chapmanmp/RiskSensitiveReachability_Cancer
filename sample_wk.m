%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Samples the disturbance at time k according to P
% INPUT:
    % ws(i) : ith possible value of wk
    % nd = length(ws), # possible values of wk
    % tick_P = [ 0, P(1), P(1)+P(2), ..., P(1)+...+P(nw-1), 1 ], nw = length(ws)
        % P(i) : probability that wk = ws(i)
 % OUTPUT:
    % wk : realization of the disturbance at time k
 % AUTHOR: Margaret Chapman
 % DATE: October 24, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wk = sample_wk( ws, nd, tick_P )

myrand = rand(1); % pseudorandom value drawn from the standard uniform distribution on (0,1)

for i = 1 : nd
    
    if myrand > tick_P(i) && myrand <= tick_P(i+1)
        
        wk = ws(i); break;
        
    end
    
end

end


