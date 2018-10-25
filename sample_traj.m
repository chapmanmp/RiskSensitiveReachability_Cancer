%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Generates one sample trajectory, cancer dynamics
% INPUT:
    % u = [u0, u1, u2] fixed control sequence, 1x3 vector
    % x = initial state, 2x1 vector
    % ak, x1s, x2s = parameters for dynamics_dt.m
    % ws, nd, tick_P = parameters to sample disturbance
    % N = 3, time horizon length
% OUTPUT: myTraj = [x0, x1, x2, x3] trajectory sample, 2x4 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function myTraj = sample_traj( u, x, ak, x1s, x2s, ws, nd, tick_P, N )

myTraj = [ x, zeros(length(x),N) ];                     % initialize trajectory

for k = 0 : N-1                                         % for each time point
    
    wk = sample_wk( ws, nd, tick_P );                   % sample the disturbance at time k according to P
    
    xk = myTraj(:,k+1);                                 % state at time k
    
    uk = u(k+1);                                        % control at time k
    
    xkPLUS1 = dynamics_dt( xk, uk, wk, ak, x1s, x2s );  % get next state realization, SNAPS TO GRID ON BDRY
    
    myTraj(:,k+2) = xkPLUS1;                            % state at time k+1
    
end

end
    
    