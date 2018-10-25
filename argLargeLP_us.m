%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Returns optimal argument to approximate for each i,j
%                   max_R\inRiskEnvelope { E[ R*J_k+1(x_k+1,ls(i)*R) | x_k, ls(i), us(j) ]
%              Uses Chow 2015 linear interpolation method on confidence level
%              Uses change of variable, Z := y*R
% INPUT:
    % f_full = [P/ls(1); P/ls(2); ...; P/ls(nl)] repeated nu times, column vec
    % bigA (matrix), bigb (col vector) : encodes linear interpolation of y*J_k+1( x_k+1, y ) vs. y given x, for all us
    % P(i): probability that w_k = ws(i)
    % ls(i): ith confidence level
    % nl = number of confidence levels
    % nd = number of disturbance values, length of P
    % nu = number of control options 
% OUTPUT:
    % |-------us(1)---------|       |--------us(nu)--------|
    % [t1,us(1);...;tnl,us(1); ... ;t1,us(nu);...;tnl,us(nu)] (col vector) that maximizes 
    % |--------------------us(1)----------------------|         |--------------------us(nu)-------------------|
    % (P/ls(1))'*t1,us(1) + ... + (P/ls(nl))'*tnl,us(1) + ... + (P/ls(1))'*t1,us(nu)+...+(P/ls(nl))'*tnl,us(nu)
        % subject to constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tStar = argLargeLP_us( f_full, bigA, bigb, P, ls, nl, nd, nu )
%big_size = max(max(abs(bigb)))/2 + max(max(abs(bigA)))/2; myfactor = max(P)/big_size; % gets constraints on similar scales
cvx_solver mosek;
for j = 1:2
    cvx_begin quiet
        if j==1, cvx_precision best; else, cvx_precision default; end

        variables Z(nd,nl*nu) t(nl*nd*nu,1)

        maximize( f_full' * t )  % (P/ls(i))' * t   
        subject to
            %(myfactor*bigA) * vec(Z) + (myfactor*bigb) >= myfactor*vec( repmat(t', nl-1, 1) );
            bigA * vec(Z) + bigb >= vec( repmat(t', nl-1, 1) );
            P' * Z == repmat( ls', 1, nu ); %ls is a column vector
            Z <= 1;
            Z >= 0;
    cvx_end

    if strcmpi(cvx_status, 'Solved') && ~isinf(cvx_optval) && ~isnan(cvx_optval)
        tStar = t; break;
    elseif j == 2
        error('maxExp.m: cvx not solved.');
    end
end

end


