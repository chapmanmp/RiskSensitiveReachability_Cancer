%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Returns all possible control sequences, N = 3
% INPUT:
    % us = possible controls at each time point
    % nu = length(us)
    % cardU = number of possible control sequences
% OUTPUT: 
    % U{i} = ith possible control sequence of the form [u0, u1, u2], 1x3 vector
% Author: Margaret Chapman
% Date: October 24, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function U = getAllControlSeq( us, nu, cardU )

U = cell(cardU,1);

s = 1;
for i = 1 : nu
    for j = 1 : nu
        for k = 1 : nu
            U{s} = [us(i) us(j) us(k)];
            s = s + 1;
        end
    end
end
    

    
    





