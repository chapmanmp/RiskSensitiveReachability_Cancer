%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Applies penalties to large # live cancer cells & large toxicity magnitude
% INPUT: State xk, or state trajectory (x0, x1, ..., xN)
    % state is [# live cells, phenotype 1; ...; # live cells, phenotype last; # dead or dying cells; toxicity]
% OUTPUT: weighted linear combination of constraint violation penalties
% AUTHOR: Margaret Chapman
% DATE: October 11, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gx = constraint_violation_penalty( x )

Tpop = 200;                                     % tolerable number of live cancer cells

xpop = x(1:end-2);                              % portion of the state with live cancer cells

Ttox = 20;                                      % tolerable toxicity magnitude

xtox = x(end);                                  % portion of the state with toxicity magnitude

gx = sum(xpop)-Tpop + Tpop/Ttox*(xtox-Ttox);    % one-sided weighted signed distance
