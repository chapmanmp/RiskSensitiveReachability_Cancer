%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Returns a proxy for toxicity of the drug applied at time k
% INPUT:
    % uk : drug applied at time k (uk = 1 is DMSO, uk = 2 is Tram, uk = 3 is BEZ, uk = 4 is Comb)
% OUTPUT: a proxy for toxicity of the drug applied at time k
% NOTE: proxy could also be the inverse of the 1-norm of the live-cell portion of A*
%           (higher toxicity for drugs that make pop shrink faster)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tox_k = getToxicityProxy( uk )

if uk == 1                  % toxicity of DMSO is 0
    
    tox_k = 0;      
    
elseif uk == 2 || uk == 3   % toxicity of Tram or BEZ is 1
    
    tox_k = 1;  
    
elseif uk == 4              % toxicity of Comb is 2
    
    tox_k = 2;  
    
else
    
    error('error in get_toxicity_proxy.m');
    
end

end