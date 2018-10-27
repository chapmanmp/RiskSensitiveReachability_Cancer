%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Applies penalties to large # live cancer cells & large toxicity magnitude
% INPUT: x = [ # live cancer cells; toxicity magnitude (proxy) ], 2x1 vector
% OUTPUT: weighted linear combination of constraint violation penalties
% AUTHOR: Margaret Chapman
% DATE: October 22, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gx = constraint_violation_penalty( x )

Tpop = 100;                                     % tolerable number of live cancer cells (small relative to 1000)

Ttox = 5;                                       % tolerable toxicity magnitude (2 applications max of combo)

%gx = x(1)-Tpop + Tpop/Ttox*(x(2)-Ttox);         % one-sided weighted signed distance

%gx = gx/1000;                                   % scale it down to accomomdate large m

% max is about 1.6

maxx1 = 999;
maxx2 = 40;

gx = mean([(x(1)-Tpop)/(maxx1-Tpop), (x(2)-Ttox)/(maxx2-Ttox)]); % max is 1

