%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Sets up cancer treatment example using PLoS model id results
    % [k, k+1) is 3 days
    % DMSO, Tram, BEZ, or Comb must be applied at each 3-day time point
    % the problem is to design an optimal schedule that reduces the risk of cancer growth and toxicity
    % the proxy for toxicity at time k is total drug concentration up to time k
% AUTHOR: Margaret Chapman
% DATE: October 10, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('model_id_K14VIMK19_july31.mat'); % PLoS model id results
% dynamics matrices for [k, k+1) duration of 12h are A_STAR{1}-DMSO, A_STAR{2}-Tram, A_STAR{3}-BEZ, A_STAR{4}-Comb

T = 6; % number of 12h intervals in 3 days

Ak = cell(N_AGENT, 1); for d = 1 : N_AGENT, Ak{d} = (A_STAR{d})^T; end 
% dynamics matrices for [k, k+1) duration of 3 days

