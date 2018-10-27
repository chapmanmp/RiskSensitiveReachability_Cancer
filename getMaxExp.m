%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: For each i,j approximates max_R\inRiskEnvelope { E[ R*J_k+1(x_k+1,ls(i)*R) | x_k, ls(i), us(j) ] }
%              Uses Chow 2015 linear interpolation method on confidence level
%              Uses change of variable, Z := y*R 
% INPUT: 
    % J_k+1: optimal cost-to-go at time k+1, matrix
    % x: state at time k, 2x1 vector
    % us: set of controls, 4x1 vector
    % nu: # controls
    % x1s: x1 grid points, row vector
    % nx1 = length(x1s)
    % x2s: x2 grid points, row vector
    % ls(j): jth discrete confidence level
    % nl: # confidence levels
    % ws(i): ith possible value of multiplicative noise wk
    % nd: # disturbance values
    % P(i): probability that w_k = ws(i)
    % ak{d}: net (noise-free) growth rate due to drug d applied at time k
        % d = 1 (DMSO), d = 2 (Tram), d = 3 (BEZ), d = 4 (Combo)
    % BIGVAL: value assigned to J_k+1(x_k+1,y) if xk+1 falls outside grid
% OUTPUT: bigexp(i,j) ~= max_R\inRiskEnvelope { E[ R*J_k+1(x_k+1,ls(i)*R) | x_k, ls(i), us(j) ] }
% AUTHOR: Margaret Chapman
% DATE: October 23, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bigexp = getMaxExp( J_kPLUS1, x, us, nu, x1s, nx1, x2s, ls, nl, ws, nd, P, ak, BIGVAL )

f_full = zeros(nd,nl); for i = 1:nl, f_full(:,i) = P/ls(i); end; f_full = vec(f_full); f_full = repmat(f_full, nu, 1);

nrows = nl*nd*(nl-1); bus = zeros(nrows,nu); Aus = []; 
for j = 1 : nu
     % encodes linear interpolation of y*J_k+1( x_k+1, y ) versus y, given us(j) and x
    [ Au, bu ] = getLMI( x, us(j), ws, nd, x1s, nx1, x2s, ls, nl, J_kPLUS1, ak, BIGVAL );

    for i = 1 : nl, Aus = blkdiag(Aus, Au); end
    
    bus(:,j) = repmat(bu, nl, 1);

end
bus = vec(bus);

%[tStar, bigexp] = argLargeLP_2( f_full, A, b, P, ls, nl, nd ); % column vector
%tStar = argLargeLP( f_full, A, b, P, ls, nl, nd ); % column vector
tStar = argLargeLP_us( f_full, Aus, bus, P, ls, nl, nd, nu ); % column vector

%bigexp = zeros(nl,1);
bigexp = zeros(nl,nu);

for j = 1 : nu
    
    tStar_j = tStar( 1 + (j-1)*nl*nd : j*nl*nd ); % extract optimal arg for us(j)

    for i = 1 : nl,  bigexp(i,j) = (P/ls(i))' * tStar_j( (i-1)*nd + 1 : i*nd ); end % extract optimal value for us(j), ls(i)
    
end

end


