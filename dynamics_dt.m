%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Defines discrete-time phenotypic-state dynamics
    % duration between k and k+1 is 3 days
% INPUT:
    % xk : [# live cancer cells, phenotype 1; ...; # live cancer cells, phenotype 4; total drug penalty] at time k
    % uk : drug applied at time k
    % wk : multiplicative error at time k
    % A_STAR{uk} : dynamics matrix for drug uk
% OUTPUT:
    % xkPLUS1 : [# live cancer cells; total drug penalty] at time k+1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xkPLUS1 = dynamics_dt( xk, uk, wk, A_STAR )

A = A_STAR{uk};

T = 6; % number of 12-hour increments in 3 days

xkPLUS1 = wk * A^T * xk; 
