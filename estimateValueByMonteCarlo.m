%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Estimates min_(u0, u1, u2) CVaR_alpha[ beta*exp(m*g(x0)) + ... + beta*exp(m*g(x3)) | x0 = x, (u0, u1, u2) ] 
%              via Monte Carlo, given state x and confidence level alpha
%              for cancer treatment example, N = 3
% INPUT: 
    % x = initial state
    % alpha = initial confidence level
    % U{i} = ith possible control sequence of the form [u0, u1, u2]
    % cardU = # possible control sequences
    % ws(i) = ith possible value of wk multiplicative modeling error
    % nd = length(ws)
    % tick_P = [ 0, P(1), P(1)+P(2), ..., P(1)+...+P(nd-1), 1 ]; P(i) = probability that wk = ws(i)
    % m, beta = parameters for stage_cost.m
    % N = 3, length of time horizon
    % ak, x1s, x2s = parameters for dynamics_dt.m
% OUTPUT: 
    % estimate of minimum CVaR_alpha cost starting from state x
    % est_uStar = a control sequence that achieved the estimated minimum         
% Author: Margaret Chapman
% Date: October 24, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ est_val, est_uStar ] = estimateValueByMonteCarlo( x, alpha, U, cardU, ws, nd, tick_P, m, beta, N, ak, x1s, x2s )

% # trials   
nt = 100000; est_cvar_u = zeros(cardU, 1);

for s = 1 : cardU % for each control sequence
    
    sample_costs = zeros(nt, 1);
    
    for q = 1 : nt
                                                                          % 2x4 matrix
        myTraj = sample_traj( U{s}, x, ak, x1s, x2s, ws, nd, tick_P, N ); % myTraj = [x0, x1, x2, x3] trajectory sample
        
        sample_costs(q) = soft_max_over_traj( myTraj, N, m, beta ) + 10^(-8)*randn(1);
        % cost of sampled trajectory + small noise to make it continuous
        % max{sample_cost} = (N+1)*beta*exp(m*1.6) = 3.6*10^4, since m = 10, beta = 10^(-3), N = 3
        % in RSReachability_Sum2018\MonteCarlo_CVaR_pond.m, max{sample_cost}=1.6*10^5 & we chose std dev = 10^(-7) 
    end
    
    var = quantile( sample_costs, 1-alpha ); % VaR_y[Z] is the (1-y)-quantile of the distribution of Z
    
    est_cvar_u(s) = estimateCVaR( sample_costs, alpha, var ); % requires continuous empirical distribution
  % est_cvar_u(s) ~= CVaR_alpha[ beta*exp(m*g(x0)) + ... + beta*exp(m*g(x3)) | x0 = x, U{s} = (u0^s, u1^s, u2^s) ]

    disp(['s = ', num2str(s)]);
    
end

[ est_val, best_s ] = min( est_cvar_u ); % choose the smallest value over all control sequences

est_uStar = U{best_s};
        
end