%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Returns value function at step N (initial value function)
% INPUT:
    % xs{i} = ith discrete state
    % nx = # discrete states
    % m, beta = stage cost parameters
    % nl = # discrete confidence levels
% OUTPUT:
    % JN(j,i) = value at time step N for state xs{i} and confidence level ls(j)
%             = beta*exp(g(xs{i}))
% Author: Margaret Chapman
% Date: October 22, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J_N = init_val_func( xs, nx, m, beta, nl )

val_states = zeros( 1, nx );

for i = 1 : nx, val_states(i) = stage_cost( xs{i}, m, beta ); end

J_N = repmat( val_states, nl, 1 ); % Initial value function, JN(x,y) = beta*exp(g(x)) for each y

