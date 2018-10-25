%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Estimates min_(u0, u1, u2) CVaR_ls(j)[ beta*exp(m*g(x0)) + ... + beta*exp(m*g(x3)) | x0 = xs{i}, (u0, u1, u2) ] 
%              via Monte Carlo, given state x and confidence level alpha
%              for cancer treatment example, N = 3
% NOTE: Compare output to Js{1}(j,i) from Main_DynProgram.m
% Author: Margaret Chapman
% Date: October 24, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clearvars; clc;

setup_example; if N ~= 3, error('N must be 3!'); end;

cardU = nu^N;                           % cardU = # possible control sequences

U = getAllControlSeq( us, nu, cardU );  % U{i} = ith possible control sequence of the form (u0, u1, u2)

i = 6100; x = xs{i};                    % fix initial state 
 
j = 8; alpha = ls(j);                   % fix initial confidence level

[ est_J0ji, est_uStar ] = estimateValueByMonteCarlo( x, alpha, U, cardU, ws, nd, tick_P, m, beta, N, ak, x1s, x2s );
% est_J0ji ~= min_(u0, u1, u2) CVaR_ls(j)[ beta*exp(m*g(x0)) + ... + beta*exp(m*g(x3)) | x0 = xs{i}, (u0, u1, u2) ] 
% est_uStar = a control sequence that achieves the estimated minimum

% Js{1}(j,i) 
