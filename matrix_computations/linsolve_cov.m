function Y = linsolve_cov(C,Y,type)
% X = linsolve_cov(C,Y)
%
% Solves C*X = Y, where C is symmetric positive definite matrix.

Y = linsolve_chol(chol(C),Y);
