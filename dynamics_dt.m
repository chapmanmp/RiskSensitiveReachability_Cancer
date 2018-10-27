%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Defines discrete-time phenotypic-state dynamics
    % duration between k and k+1 is 2 days
    % state is [# live cancer cells; toxicity magnitude (proxy) ]
% INPUT:
    % xk : state at time k, 2x1 vector
    % uk : drug applied at time k, scalar
    % wk : multiplicative modeling error on [k, k+1), scalar
    % ak{uk} : net growth rate for 2-day interval, drug uk
    % x1s: x1 grid points, row vector
    % x2s: x2 grid points, row vector
% OUTPUT:
    % xkPLUS1 : state at time k+1, 2x1 vector
% AUTHOR: Margaret Chapman
% DATE: October 23, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xkPLUS1 = dynamics_dt( xk, uk, wk, ak, x1s, x2s )

xkPLUS1 = xk;             

xkPLUS1(1) = wk * ak{uk} * xk(1);  
% number of live cancer cells at time k+1

xkPLUS1(2) = xk(2) + getToxicityProxy(uk); 
% toxicity magnitude (proxy) at time k+1 due to drug uk applied at time k (and previous drugs)

% snap to fixed point outside of grid - over long time horizons cost may become too large for numerical stability
f = 2;
if xkPLUS1(1) > f*max(x1s), xkPLUS1(1) = f*max(x1s); end

if xkPLUS1(2) > f*max(x2s), xkPLUS1(2) = f*max(x2s); end

end





