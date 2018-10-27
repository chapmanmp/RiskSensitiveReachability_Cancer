%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Computes the parameters for linear matrix inequality given x and u
% INPUT: 
    % x: state at time k, 2x1 vector
    % u: control at time k, scalar
    % ws(i): ith realization of disturbance, scalar
    % nd: # disturbance values
    % x1s: x1 grid points, row vector
    % nx1 = length(x1s)
    % x2s: x2 grid points, row vector
    % ls(j): jth discrete confidence level, scalar
    % nl: # confidence levels
    % J_k+1(j,i): cost-to-go for sub-problem that starts at time k+1, state xs{i}, confidence level ls(j)
    % ak{d}: net (noise-free) growth rate due to drug d applied at time k
        % d = 1 (DMSO), d = 2 (Tram), d = 3 (BEZ), d = 4 (Combo)
    % BIGVAL: value assigned to J_k+1(x_k+1,y) if xk+1 falls outside grid
% OUTPUT: A = blkdiag( A1, ..., And ), b = [b1; ...; bnd]
%             |A1 0  .. 0   |                |b1 |
%           = |0  ..... 0   |              = |...|
%             |0  0  .. And |                |bnd|
% NOTE:
    % Ai & bi are column vectors that encode the linear interpolation of y*J_k+1( x_k+1, y ) vs. y, given x and u
        % at the ith realization of x_k+1 = dynamics_dt( x, u, ws(i), ak )
    % max_t,y { t | A1(j)*y + b1(j) >= t, confidence level line segment j } is equivalent to 
        % max_y { g(y) := min_j A1(j)*y + b1(j), confidence level line segment j }                                          
    % g(y) = linear interpolation of y*J_k+1(x,y) vs. y, at fixed x (and u); concave & piecewise linear in y
    % uses Chow, et al. NIPS 2015 to manage continuous confidence level
    % uses bilinear interpolation to manage continuous state space
% AUTHOR: Margaret Chapman
% DATE: October 23, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ A, b ] = getLMI( x, u, ws, nd, x1s, nx1, x2s, ls, nl, J_kPLUS1, ak, BIGVAL )

A = []; b_mat = zeros(nl-1,nd); % to contain [b1 b2 ... bnd]

for i = 1 : nd                  % for each disturbance realization

    x_kPLUS1 = dynamics_dt( x, u, ws(i), ak, x1s, x2s ); % get next state realization, 2x1 vector (SNAPS TO GRID ON BDRY)
    
    Ai = zeros(nl-1,1); bi = zeros(nl-1,1);              % one entry per confidence level line segment
        
    for j = nl-1: -1: 1 % for each confidence level line segment, [l_j+1, l_j], e.g., ls = [l_1 = 0.95, l_2 = 1/2, l_3 = 0.05] 
                        % [l_3, l_2] = [0.05, 1/2] 
                        % [l_2, l_1] = [1/2, 0.95]
                          
        J_jPLUS1 = interp2( x1s, x2s, vec2mat(J_kPLUS1(j+1,:),nx1), x_kPLUS1(1), x_kPLUS1(2), '*linear', BIGVAL ); 
        % approximates J_k+1(x_k+1, ls(j+1)) using bilinear interpolation via J_k+1(x_g, ls(j+1)), where x_g \in grid
        % x1s, x2s are (and must be) evenly spaced and monotonic
        % assigns value BIGVAL if x_kPLUS1(i) is outside of xis
              
        J_j = interp2( x1s, x2s, vec2mat(J_kPLUS1(j,:),nx1), x_kPLUS1(1), x_kPLUS1(2), '*linear', BIGVAL );
        % approximates J_k+1(x_k+1, ls(j))
        
        Ai(j) = ( ls(j)*J_j - ls(j+1)*J_jPLUS1 )/( ls(j)-ls(j+1) ); 
        % approx. slope of jth line segment of linear_interp( y*J_k+1( x_k+1, y ) vs. y ) 
        
        bi(j) = ls(j+1) * (J_jPLUS1 - Ai(j));                       
        % approx. y-int of jth line segment of linear_interp( y*J_k+1( x_k+1, y ) vs. y )
    end
    
    A = blkdiag( A, Ai ); b_mat(:,i) = bi; %fillup matrix column by column   
    
end

b = vec(b_mat); % vectorizes matrix

end
