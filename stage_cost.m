%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Defines the stage cost as an exponential of the constraint violation penalty
% INPUT: 
    % x: 2x1 state vector
    % m: greater than or equal to 1, soft-max parameter
    % beta: scaling parameter so LP CVX solver works for large m
% OUTPUT: Stage cost
% AUTHOR: Margaret Chapman
% DATE: October 22, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = stage_cost( x, m, beta ) 

gx = constraint_violation_penalty( x );  % constraint violation penalty, vector or matrix

c = beta * exp( m*gx );