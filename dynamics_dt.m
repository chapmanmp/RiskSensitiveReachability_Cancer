%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Defines discrete-time phenotypic-state dynamics
    % duration between k and k+1 is 3 days
    % state is [# live cells, phenotype 1; ...; # live cells, phenotype last; # dead or dying cells; toxicity]
% INPUT:
    % xk : state at time k
    % uk : drug applied at time k
    % wk : multiplicative modeling error on [k, k+1)
    % Ak{uk} : dynamics matrix for 3-day duration, drug uk 
% OUTPUT:
    % xkPLUS1 : state at time k+1
% AUTHOR: Margaret Chapman
% DATE: October 10, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xkPLUS1 = dynamics_dt( xk, uk, wk, Ak )

xkPLUS1 = zeros(size(xk));             

xkPLUS1(1:end-1) = wk * Ak{uk} * xk(1:end-1);  % cell population portion of the state

xkPLUS1(end) = xk(end) + getToxicityProxy(uk); % toxicity at time k+1 due to drug uk applied at time k (and previous drugs)

