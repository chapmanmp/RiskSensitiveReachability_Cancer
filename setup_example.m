%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Sets up cancer treatment example using PLoS model id results
    % [k, k+1) is 2 days
    % DMSO, Tram, BEZ, or Comb must be applied at each 2-day time point
    % the problem is to design an optimal schedule that reduces the risks of cancer growth and toxicity
    % the proxy for toxicity at time k is total drug concentration up to time k
% AUTHOR: Margaret Chapman
% DATE: October 22, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('model_id_K14only_aug6.mat');      % PLoS model id results

clearvars -except A_STAR N_AGENT        % dynamics matrices for 12h duration
                                        % A_STAR{1}(DMSO), A_STAR{2}(Tram), A_STAR{3}(BEZ), and A_STAR{4}(Comb)

T = 4;                                  % number of 12h intervals in 2 days

ak = cell(N_AGENT, 1); 
for d = 1 : N_AGENT
                                        % net growth rate for 12h interval
    mu_d = norm( A_STAR{d}(1:end-1, 1:end-1), 1 );  
    
    ak{d} = mu_d^T;                     % net growth rate for 2-day interval

end 

x1s = 0 : 200 : 1000; nx1 = length(x1s); % Discretized x1's    [ # live cancer cells in population ]

x2s = 0 : 10 : 40; nx2 = length(x2s);   % Discretized x2's    [ toxicity magnitude (proxy) ]

nx = nx1*nx2;                          % Number discrete states in total

xs = cell(nx,1); xs_array = cell(nx1,nx2); n = 1;

for j = 1 : nx2
    for i = 1 : nx1
        xs{n} = [x1s(i); x2s(j)];
        xs_array{i,j} = xs{n};
        n = n + 1;
    end 
end

ls = (0.95 : -0.1 : 0.05)';             % Discretized confidence levels
nl = length(ls);

[ X1, L ] = meshgrid( x1s, ls' );       % for plotting at fixed x2

N = 10;                                 % Time horizon: {0, 1, 2, ..., N} = {0, 2, 4, ..., 20} days

ws = [0.5, 1, 2];                       % ws(i): ith possible value of wk multiplicative modeling error

P = getProbDist( ws, 1.5 );             % P(i): probability that wk = ws(i), choose mean = 1.5 to be de-stabilizing

m = 10;                                 % soft-max parameter

beta = 10^(-3);                         % scaling parameter 

