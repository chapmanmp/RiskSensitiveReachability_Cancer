%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Performs the CVaR Bellman backwards recursion, uk \in {0,1}
% INPUT: 
    % J_k+1: optimal cost-to-go at time k+1, matrix
    % xs{i}: ith discrete state
    % nx: number of discrete states
    % x1s: x1 grid points, row vector
    % nx1 = length(x1s)
    % x2s: x2 grid points, row vector
    % ls(j): jth discrete confidence level
    % nl: # confidence levels
    % ws(i): ith possible value of multiplicative noise wk
    % nd: # disturbance values
    % us: set of controls; 1: DMSO, 2: Tram, 3: BEZ, 4: Combo
    % nu: # controls
    % P(i): probability that wk = ws(i)
    % m, beta: parameters for stage_cost.m
    % ak{d}: net (noise-free) growth rate due to drug d applied at time k
        % d = 1 (DMSO), d = 2 (Tram), d = 3 (BEZ), d = 4 (Combo)
% OUTPUT: 
    % J_k(j,i): approx. best cost-to-go for sub-problem that starts at time k, state xs{i}, confidence level ls(j)
    % mu_k(j,i): approx. best drug input at time k, state xs{i}, confidence level ls(j)
% AUTHOR: Margaret Chapman
% DATE: October 23, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ J_k, mu_k ] = CVaR_Bellman_Backup( J_kPLUS1, xs, nx, x1s, nx1, x2s, ls, nl, ws, nd, us, nu, P, m, beta, ak )

J_k = J_kPLUS1; mu_k = J_kPLUS1; % initialization

for i = 1 : nx      % <--x's change along columns of J_k-->
    
    x = xs{i};
    
    maxExp_us = getMaxExp( J_kPLUS1, x, us, nu, x1s, nx1, x2s, ls, nl, ws, nd, P, ak ); 
    %maxExp_us(i,j) = maxExp, given state x, confidence level ls(i), and control us(j)
    
    %maxExp_u1 = maxExp_pond( J_kPLUS1, x, us(1), xs, ls, ws, P, dt, area_pond ); 
    %maxExp_u2 = maxExp_pond( J_kPLUS1, x, us(2), xs, ls, ws, P, dt, area_pond );
    %[ optExp, optInd ] = min( [maxExp_u1, maxExp_u2], [], 2 ); % col vector, 1 entry per confidence level
   
    [ optExp, optInd ] = min( maxExp_us, [], 2 ); % col vector, 1 entry per confidence level
    
    J_k(:,i) = stage_cost( x, m, beta ) + optExp;
    
    mu_k(:,i) = us(optInd); % need to check
    
end

end