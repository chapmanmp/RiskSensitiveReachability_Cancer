%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Returns sum of exponentials, beta*exp(m*g(x0)) + beta*exp(m*g(x1)) + ... + beta*exp(m*g(xN))
% INPUT:
    % myTraj = [x0, x1, ..., xN] trajectory sample, 2x(N+1) matrix
    % N = time horizon length
    % m, beta = parameters for stage_cost.m
% OUTPUT: weighted sum of exponentials over a given trajectory
% Author: Margaret Chapman
% Date: October 24, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sample_cost = soft_max_over_traj( myTraj, N, m, beta )

sample_cost = 0;

for k = 0 : N
    
    xk = myTraj(:,k+1); % state at time k
    
    sample_cost = sample_cost + stage_cost( xk, m, beta ); 
    
end