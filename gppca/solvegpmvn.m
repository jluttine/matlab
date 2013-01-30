
function [sqrtCov, Cov, mu] = solvegpmvn(K, A, B)
% sqrtCov = solvegpmvn(K, A)
% [sqrtCov, Cov] = solvegpmvn(K, A)
% [sqrtCov, Cov, mu] = solvegpmvn(K, A, B)
%
% Solves the set of equations
%  Cov^(-1) = K^(-1) * A * K^(-1)
%  mu = Cov * K^(-1) * B
% or equivalently for Cov,
%  mu = K * A^(-1) * B
%  Cov = K * A^(-1) * K
% using Cholesky decomposition. K and A should be symmetric positive
% definite matrices of the same size. B should be a column vector.
%
% sqrtCov is such that Cov = sqrtCov'*sqrtCov

U = safechol(A); % A=U'*U, i.e., A^(-1) = U^(-1) * U^(-1)'
opts.UT = true;
opts.TRANSA = true;
sqrtCov = linsolve(U,K,opts); % sqrtCov = U^(-1)' * K
if nargout >= 2
  Cov = sqrtCov'*sqrtCov;
end
if nargout >= 3
  mu = sqrtCov'*linsolve(U,B,opts); % mu = K * U^(-1) * U^(-1)' * B
end
